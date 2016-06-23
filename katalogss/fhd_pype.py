import numpy as np
import os
from scipy.io import readsav
from astropy.io import fits

global fhd_base 
fhd_base = '/nfs/eor-09/r1/djc/EoR2013/Aug23/'

def set_fhd_base(path):
    global fhd_base
    if not path.endswith('/'): path+='/'
    fhd_base = path
    return


def get_obslist(fhd_run):
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base,fhd_run)
    obs=os.listdir(decon_dir)
    obs=[o[:10] for o in obs if o.endswith('fhd.sav')]
    obs.sort() 
    return obs


def read_sourcelist(fhdsav, tag='component_array', keys=('id','x','y','ra','dec','flux','gain','alpha','freq')):
	cat = readsav(fhdsav)[tag]
        items = [cat.id,cat.x,cat.y,cat.ra,cat.dec,np.vstack(cat.flux).T['i'][0],cat.gain,cat.alpha, cat.freq]
        items = [item.astype(np.float64) for item in items]
        cat = dict(zip(keys,items))
	return cat


def gen_cal_cat(cat, freq=180., alpha=-0.8, file_path='catalog.sav'):
    cat_keys=cat.keys()
    cat['id']=np.argsort(cat['flux'])[::-1]
    nsrcs = len(cat['id'])
    if 'alpha' in cat_keys:pass
    else: cat['alpha'] = alpha*np.ones(nsrcs)
    if 'freq' in cat_keys:pass
    else: cat['freq']=float(freq)*np.ones(nsrcs)
    import pidly
    idl=pidly.IDL()
    idl.id=cat['id']
    idl.ra=cat['ra']
    idl.dec=cat['dec']
    idl.n=len(cat['id'])
    idl.freq=float(freq)
    idl.alpha=cat['alpha']
    idl.fluxi=cat['flux']
    idl('sl=source_comp_init(n_sources=n,ra=ra,dec=dec,freq=freq,flux=fluxi,alpha=alpha,id=id)')
    idl('generate_calibration_catalog,sl,file_path=\'%s\''%file_path)
    idl.close()
    return cat


def fetch_comps(fhd_run, obsids=None, cache=True):
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base,fhd_run)
    meta_dir='%sfhd_%s/metadata/'%(fhd_base,fhd_run)
    if obsids is None: obsids = get_obslist(decon_dir)
    comps={} 
    for o in obsids:
        fhdsav=decon_dir+o+'_fhd.sav'
	print 'fetching data for obsid: %s'%o
	comps[o] = read_sourcelist(fhdsav)
    return comps


def fetch_meta(fhd_run, obsids=None):
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base,fhd_run)
    meta_dir='%sfhd_%s/metadata/'%(fhd_base,fhd_run)
    if obsids is None: obsids = fp.get_obslist(decon_dir)
    meta = {}
    for o in obsids:
        params = readsav(decon_dir+o+'_fhd_params.sav')['fhd_params']       
        metaobs = readsav('%s%s_obs.sav'%(meta_dir,o))['obs']
        meta[o] = {'n_iter':params.n_iter[0],'det_thresh':params.detection_threshold[0],'beam_thresh':params.beam_threshold[0],'max_bl':metaobs.max_baseline[0],'freq':metaobs.freq_center[0],'degpix':metaobs.degpix[0]}
        meta[o]['beam_width'] = meta[o]['max_bl']**-1 * 180./np.pi
    return meta


def pixarea_maps(fhd_run, obsids=None, map_dir='area_maps/', recalculate=False):
    if not os.path.exists(map_dir): os.system('mkdir %s'%map_dir)
    if not map_dir.endswith('/') map_dir += '/'
    if obsids is None: obsids = fp.get_obslist(decon_dir)
    calcobs = [o for o in obsids if recalculate or not os.path.exists(map_dir+o+'_amap.fits')]
    if len(calcobs)>0:
        import pidly
        idl=pidly.IDL()
        idl.fhddir='%sfhd_%s/'%(fhd_base,fhd_run)
        idl.mapdir = map_dir
        for o in calcobs:
            idl.obsid=o
            commands = ['restore,fhddir+\'metadata/\'+obsid+\'_obs.sav\'',
                    'area_map=pixel_area(obs)',
                    'beam_width = 1/obs.max_baseline',
                    'beam_area = !pi*beam_width^2/(4*alog(2))',
                    'area_map /= beam_area',
                    'writefits, mapdir+obsid+\'_amap.fits\',area_map']
            for c in commands:idl(c)
        idl.close()
    amaps=dict([(o,fits.open(map_dir+o+'_amap.fits')[0].data) for i,o in enumerate(obsids)])
    return amaps


def get_maps(fhd_run, obsids, imtype):
    map_dir='%sfhd_%s/output_data/'%(fhd_base,fhd_run)
    imgs = {}
    for o in obsids:
        fname = map_dir+o+'_'+imtype+'.fits'
        if os.path.exists(fname): 
            hdu = fits.open(fname)[0]
            imgs[o] = hdu
        else: imgs[o] = None
    return imgs

