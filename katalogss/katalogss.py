'''
Generate source catalogs from FHD deconvolution output.
'''


import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.modeling import models, fitting
from sklearn.cluster import DBSCAN
from kg_utils import centroid, sigma_clip


def clip_comps(comps, nmax=None):
     ''' 
     Removes flagged components.

     Parameters
     ----------
     comps: `pandas.DataFrame`
         Deconvolved souce components. A \'flag\' column is expected.
     nmax: `int`
         Maximum number of comonents to keep.
     '''

     comps = comps[comps.flag==0]
     if nmax is not None: comps = comps.iloc[:nmax]
     return comps

#def optimize_eps():


def cluster_sources(sources, eps, min_samples=1):
     ''' 
     Clusters components into spatially isolated sources. 

     Parameters
     ----------
     sources: `pandas.DataFrame`
         Sources or components to be clustered. \'ra\' and \'dec\' columns are expected.
     eps: `float`
         Radius within which two points are considered part of the same source. 
     min_samples: `int`, optional
         Minimum number of sources needed to define a cluster. The default is 1.
     '''

     a,b = sources[['ra','dec']].values.T
     X= np.array([a,b]).T
     db = DBSCAN(eps=eps, min_samples=min_samples)
     db.fit(X)
     sources.index = db.labels_
     return sources
     

def _calc_source_params(comps, labelset=None):
     ''' 
     Called by `catalog_sources` to calculate source parameters. This is currently the bottle neck function call and could use improvements. 

     Parameters
     ----------
     comps: `pandas.DataFrame`
          Components data frame indexed by cluster label. \'ra\', \'dec\', \'flux\', and \'gain\' columns are expected. 
     '''

     labelset = set(comps.index)
     keys=['id','ncomp','ra','sig_ra','dec','sig_dec','flux']
     params=[]
     # This is a major bottleneck. It takes about 30 minutes for 20000 sources.
     for i in labelset:
          ra, dec, flux, gain = comps.loc[i][['ra', 'dec', 'flux', 'gain']].values.T
          ncomp = (comps.index == i).sum()
          sra, ssra = centroid(ra, flux)
          sdec, ssdec =centroid(dec,flux)
          g = 1.-np.product(1-gain)
          sflux = np.sum(flux)/g     
          params.append([i,ncomp,sra,ssra,sdec,ssdec,sflux]) 
     return pd.DataFrame(columns=keys, data=params)


def _calc_local_rms(residual, coords, width=20, sig_clip=None):   
     '''
     Calculate the rms of residual image within a specified box region.

     Parameters
     ----------
     residual: `array`
          2D residual image array.
     coords: `tuple`
          (y,x) coordinates at which to estimate the residual image rms.
     width : `float`, optional
          Width of the box region within which to calculate the residual rms. default is 20 pixels.
     sig_clip: `float`, optional
          Rhe significance level at which to clip outliers. Default is `None`.
     '''
     
     y,x = coords
     r,mr = width/2,width%2  # The mod is added to maintain odd widths.
     res = resdiual[y-r,y+r+mr, x-r, x+r+mr].flatten()
     if sig_clip is not None: res = sigma_clip(res,sig_clip)
     rms = np.sqrt(np.mean(res**2))
     return rms


def _rms_vs_beam(residual, beam, beam_thresh=0.1): 
     '''
     Estimate the rms of residual image as a function of beam response.

     Parameters
     ----------
     residual: `array`
          2D residual image array.
     beam: `array`
          2D beam image array.
     beam_thresh : `float`, optional
          Minimum inclusive beam value.
     '''
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


def catalog_sources(comps, meta, residual, beam):
     '''
     Generate catalog of source cnadidates from components.

     Parameters
     ----------
     comps : `pandas.DataFrame`
         Deconvolved souce components. A \'flag\' column is expected.
     meta: `dict`
         FHD meta data.         
     residual: `astropy.fits.HDU` object
         Residual image in Jy/beam.
     beam: `astropy.fits.HDU` object
         Beam image array.
     '''
     cat = _calc_source_params(comps)
     wcs = WCS(beam.header)
     x,y = wcs.wcs_world2pix(cat.ra.values,cat.dec.values,1)
     cat['x'],cat['y'] = x.astype(int), y.astype(int)
     cat['beam'] = beam.data[[cat.y,cat.x]]
     cat['sig_flux'] = _rms_vs_beam(residual.data,beam.data)(cat.beam)
     cat['ext'] = (1.665 * np.sqrt(cat.sig_ra**2 + cat.sig_dec**2)) > meta['beam_width']
     return cat


