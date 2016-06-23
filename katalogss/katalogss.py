import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.modeling import models, fitting
from sklearn.cluster import DBSCAN


def clip_comps(comps):
     try: comps = comps.iloc[: np.where(comps.id==-1)[0][0]]
     except IndexError: pass
     return comps


#def optimize_eps():


def cluster_sources(sources, eps, min_samples, coords='sky'):
     if coords=='pix': a,b = sources[['x','y']].values.T
     else: a,b = sources[['ra','dec']].values.T
     a[a<180.]+=360.
     X= np.array([a,b]).T
     db = DBSCAN(eps=eps, min_samples=min_samples)
     db.fit(X)
     sources.index = db.labels_
     sources.clustered = True
     return sources


def _flux_mean(a, flux):
     mu = np.sum(a*flux)/np.sum(flux)
     sd = np.sqrt(np.sum(flux * (a-mu)**2)/np.sum(flux))
     return  mu,sd


def _calc_source_params(comps, labelset=None):
     if labelset is None: labelset = set(comps.index)
     keys=['id','ncomp','ra','sig_ra','dec','sig_dec','flux']
     params=[]
     for i in labelset: 
          ra,dec,flux,gain=comps.loc[i][['ra','dec','flux','gain']].values.T
          ncomp = (comps.index==i).sum()
          sra, ssra = _flux_mean(ra,flux)
          sdec, ssdec = _flux_mean(dec,flux)
          g = 1.-np.product(1-gain)
          sflux = np.sum(flux)/g     
          params.append([i,ncomp,sra,ssra,sdec,ssdec,sflux]) 
     return pd.DataFrame(columns=keys, data=params)


def _calc_local_rms(residual, coords, r=10, sig_clip=1):          
     y,x = coords
     res = resdiual[y-r,y+r+1, x-r, x+r+1].flatten()
     if sig_clip: res = res[np.abs(res-np.mean(res))<(sig_clip*np.std(res))]
     rms = np.sqrt(np.mean(res**2))
     return rms


def _rms_vs_beam(residual, beam, beam_thresh=0.1):  
     bins = np.arange(beam_thresh,1.01,.01)
     bin_centers = bins[:-1]+np.diff(bins)
     rms=np.zeros_like(bin_centers)
     for i,(b1,b2) in enumerate(zip(bins[:-1],bins[1:])):
          res = residual[(beam>b1)&(beam<=b2)].flatten()
          rms[i] = np.sqrt(np.mean(res**2))
     mask = ~np.isnan(rms)  
     p = models.PowerLaw1D(amplitude=np.mean(rms[mask]), x_0=1., alpha=1.) 
     p.x_0.fixed=True
     pfit = fitting.LevMarLSQFitter()
     p = pfit(p,bin_centers[mask],rms[mask])                
     return p


def catalog_sources(comps,  residual, beam): 
     cat = _calc_source_params(comps)
     wcs = WCS(beam.header)
     x,y = wcs.wcs_world2pix(cat.ra,cat.dec,1)
     cat['x'],cat['y'] = x.astype(int), y.astype(int)
     cat['beam'] = beam.data[[cat.y,cat.x]]
     cat['sig_flux'] = _rms_vs_beam(residual.data,beam.data)(cat.beam)
     return cat


