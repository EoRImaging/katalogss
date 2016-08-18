'''
Read and manipulate FHD output in Python.
'''


import os
import numpy as np
from scipy.io import readsav
from astropy.io import fits


global _fhd_base
_fhd_base = '/nfs/eor-09/r1/djc/EoR2013/Aug23/'


def fhd_base():
    'Return path to FHD output directory.'
    global _fhd_base
    return _fhd_base


def set_fhd_base(path):
    '''
    Set the FHD output directory returned by `fhd_base` method.

    Parameters
    ----------
    path: string
        The path to the FHD output directory
    '''
    global _fhd_base
    if not path.endswith('/'): path+='/'
    _fhd_base = path
    if not os.path.exists(_fhd_base): raise ValueError,'%s does not exist. Use `set_fhd_base` to update.'%_fhd_base
    return


def get_obslist(fhd_run):
    '''
    Get the list of obsids with deconvolution ouput.

    Parameters
    ----------
    fhd_run: string
        The name identifier of the FHD run, e.g. \'pac_decon_eor1_June2016\'.
    '''
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base(),fhd_run)
    obs=os.listdir(decon_dir)
    obs=[o[:10] for o in obs if o.endswith('fhd.sav')]
    obs.sort() 
    return obs


def read_sourcelist(fhdsav, tag='component_array'):
    '''
    Get the list of obsids with deconvolution ouput.

    Parameters
    ----------
    fhdsav: string    fhd_run: string

        Full path to the IDL save file containing a source list structure.
    tag: string, optional
        The tag name of the source list in the IDL structure. Defaults to \'component_array\'.
    '''
    cat = readsav(fhdsav)[tag]
    items = [cat.id, cat.x, cat.y, cat.ra, cat.dec, 
             np.vstack(cat.flux).T['i'][0],
             cat.gain,cat.alpha, cat.freq, cat.flag]
    items = [item.astype(np.float64) for item in items]
    cat = dict(zip(['id','x','y','ra','dec','flux','gain','alpha','freq',
                    'flag'],items))
    return cat


def gen_cal_cat(cat, freq=180., alpha=-0.8, file_path='catalog.sav'):
    '''
    Generate IDL structure and save file from `katalogss` catalog dict for input to FHD.

    Parameters
    ----------
    cat: dict
        The source catalog. Keys (\'ra\', \'dec\', \'flux\') are required. Keys (\'alpha\', \'freq\') are optional.
    freq: float, optional
        Frequency (MHz) assigned only if \'freq\' is not in `cat`. Defaults to 180.
    alpha: float, optional
        Spectral index assigned only if \'alpha\' is not in `cat`. Defaults to -0.8.
    file_path: string, optional
        File path passed to generate_calibration_catalog.pro. Defaults to \'catalog.sav\'.
    '''
    cat_keys=cat.keys()
    cat['id']=np.argsort(cat['flux'])[::-1]
    nsrcs = len(cat['id'])
    if 'alpha' in cat_keys:pass
    else: cat['alpha'] = alpha*np.ones(nsrcs)
    if 'freq' in cat_keys:pass
    else: cat['freq']=float(freq)*np.ones(nsrcs)
    import pidly
    idl=pidly.IDL()
    idl('!Quiet=1')
    idl.id=cat['id']
    idl.ra=cat['ra']
    idl.dec=cat['dec']
    idl.n=len(cat['id'])
    idl.freq=cat['freq']
    idl.alpha=cat['alpha']
    idl.fluxi=cat['flux']
    idl('sl=source_comp_init(n_sources=n,ra=ra,dec=dec,freq=freq,flux=fluxi,alpha=alpha,id=id)')
    idl('generate_calibration_catalog,sl,file_path=\'%s\''%file_path)
    idl.close()
    return cat


def fetch_comps(fhd_run, obsids=None):
    '''
    Return the FHD deconvolved source components.

    Parameters
    ----------
    fhd_run: string
        The name identifier of the FHD run, e.g. \'pac_decon_eor1_June2016\'.
    obsids: list-like, optional
        Obsids (as strings) to fetch data from. Defaults to all deconvolved.
    '''
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base(),fhd_run)
    meta_dir='%sfhd_%s/metadata/'%(fhd_base(),fhd_run)
    if obsids is None: obsids = get_obslist(decon_dir)
    comps={} 
    for o in obsids:
        print 'Fetching compoonent array from %s_fhd.sav'%o
        fhdsav=decon_dir+o+'_fhd.sav'
	comps[o] = read_sourcelist(fhdsav)
    return comps


def fetch_meta(fhd_run, obsids=None):
    '''
    Return meta data needed for the FHD deconvolved source components.

    Parameters
    ----------
    fhd_run: string
        The name identifier of the FHD run, e.g. \'pac_decon_eor1_June2016\'.
    obsids: list-like, optional
        Obsids (as strings) to fetch data from. Defaults to all deconvolved.
    '''
    decon_dir='%sfhd_%s/deconvolution/'%(fhd_base(),fhd_run)
    meta_dir='%sfhd_%s/metadata/'%(fhd_base(),fhd_run)
    if obsids is None: obsids = fp.get_obslist(decon_dir)
    meta = {'clustered':False}
    for o in obsids:
        params = readsav(decon_dir+o+'_fhd_params.sav')['fhd_params']       
        metaobs = readsav('%s%s_obs.sav'%(meta_dir,o))['obs']
        meta[o] = {'n_iter':params.n_iter[0],'det_thresh':params.detection_threshold[0],'beam_thresh':params.beam_threshold[0],'max_bl':metaobs.max_baseline[0],'freq':metaobs.freq_center[0],'degpix':metaobs.degpix[0]}
        meta[o]['beam_width'] = meta[o]['max_bl']**-1 * 180./np.pi
    return meta


def pixarea_maps(fhd_run, obsids=None, map_dir='area_maps/', recalculate=False):
    '''
    Return the pixel area maps and cache locally. 

    Parameters
    ----------
    fhd_run: string
        The name identifier of the FHD run, e.g. \'pac_decon_eor1_June2016\'.
    obsids: list-like, optional
        Obsids (as strings) to fetch data from. Defaults to all deconvolved.
    map_dir: string, optional
        The directory in which to cache the area maps. Defaults to \'area_maps/\'.
    recalculate: bool, optional
        If `True` and `map_dir` exists, re-run the IDL code and re-cache. 
    '''
    if not os.path.exists(map_dir): os.system('mkdir %s'%map_dir)
    if not map_dir.endswith('/'): map_dir += '/'
    if obsids is None: obsids = fp.get_obslist(decon_dir)
    calcobs = [o for o in obsids if recalculate or not os.path.exists(map_dir+o+'_amap.fits')]
    if len(calcobs)>0:
        import pidly
        idl=pidly.IDL()
        idl.fhddir='%sfhd_%s/'%(fhd_base(),fhd_run)
        idl.mapdir = map_dir
        for o in calcobs:
            idl.obsid=o
            commands = ['!quiet=1','!except=0','restore,fhddir+\'metadata/\'+obsid+\'_obs.sav\'',
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
    '''
    Return a dictionary of image maps for given obsids. Images are `astropy.fits.HDU` objects with `header` and `data` attributes.

    Parameters
    ---------
    fhd_run: string
        The name identifier of the FHD run, e.g. \'pac_decon_eor1_June2016\'.
    obsids: list-like
        Obsids (as strings) to data from.
   imtype: string
       Specifies the image type as found in the fits file name, e.g. \'uniform_Residual_I\'.
    '''
    map_dir='%sfhd_%s/output_data/'%(fhd_base(),fhd_run)
    imgs = {}
    for o in obsids:
        fname = map_dir+o+'_'+imtype+'.fits'
        if os.path.exists(fname): 
            hdu = fits.open(fname)[0]
            imgs[o] = hdu
        else: imgs[o] = None
    return imgs


