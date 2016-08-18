'''
This is a high level script designed to be run from the command line to generate source catalogs from FHD deconvolution output. It defaults to run on the MIT cluster, in which case the only required input is the string identifier for the run (e.g. \'pac_decon_eor1_June2016\'). Use the --help option to see usage.
'''


import os
import pickle
import multiprocessing as mp
from optparse import OptionParser
import warnings
import numpy as np
import pandas as pd
import fhd_pype as fp
import katalogss as kg
from kg_utils import approx_stokes_i
warnings.filterwarnings('ignore')


def get_images():
    '''
    Return the Beam_XX, Beam_YY, and Residual_I images.
    '''
    beamxx, beamyy, residual = [fp.get_maps(fhd_run, obsids=obsids, imtype=imtype) for imtype in ('Beam_XX','Beam_YY','uniform_Residual_I')]
    pix2beam = fp.pixarea_maps(fhd_run, obsids=obsids, map_dir=kgs_out+'area_maps/')
    residual = [r.data * pix2beam[o] for r,o in zip(residual,obsids)]
    return beamxx,beamyy,residual


def remove_badobs():
    '''
    Remove obsids if beam or residual data is missing.
    '''
    badobs=[o for o in obsids if None in (residual[o], beamxx[o], beamyy[o])]
    for o in badobs: obsids.remove(o)
    return obsids
    

def get_component_data():
    '''
    Return deconvolved source components.
    '''
    pool = mp.Pool()
    prcs = [pool.apply_async(fp.fetch_comps, args=(fhd_run,), kwds={'obsids': [obsid]}) for obsid in obsids]
    comps= dict([(obsid,prc.get().values()[0]) for obsid,prc in zip(obsids,prcs)])
    pool.terminate()
    meta = fp.fetch_meta(fhd_run,obsids=obsids)
    print 'Saving component data to \n%s'%comp_file
    pickle.dump([comps,meta], open(comp_file,'w'))
    return comps,meta


def _cluster_one(obsid):
    '''
    Cluster components into sources for a single obsid.
    '''
    cmps = pd.DataFrame(comps[obsid])
    cmps = kg.clip_comps(cmps,nmax=max_comps)
    cmps.ra[cmps.ra<shift_ra]+=360. 
    print 'Clustering obsid %s...'%obsid
    cmps = kg.cluster_sources(cmps,  eps_factor * meta[obsid]['beam_width'])
    cmps[cmps.ra>=360.]-=360.
    return cmps


def _catalog_one(obsid):
    '''
    Generate catalog from clustered components for a single obsid.
    '''
    beam = approx_stokes_i(beamxx[obsid], beamyy[obsid])
    print 'Cataloging obsid %s...'%obsid
    catalog = kg.catalog_sources(comps[obsid], meta[obsid], residual[obsid], beam)
    return catalog


def cluster_components():
    '''
    Batch cluster components to sources for all obsids.
    '''
    pool = mp.Pool()
    comps =  dict(zip(obsids,pool.map(_cluster_one,obsids)))
    pool.terminate()
    meta['clustered']=True
    meta['eps_factor']=eps_factor
    print 'Saving clustered component data to \n%s'%comp_file
    pickle.dump([comps,meta],open(comp_file,'w'))
    return comps,meta


def generate_catalog():
    '''
    Batch catalog sources for all obsids.
    '''
    pool = mp.Pool()
    cats =  dict(zip(obsids,pool.map(_catalog_one,obsids)))
    pool.terminate()
    print 'Saving catalogs to \n%s'%cat_file
    pickle.dump(cats,open(cat_file,'w'))
    return


def get_args():
    '''
    Parse command line arguments. Note optparse is deprecated. This should be re-written with argparse.
    '''
    parser = OptionParser(usage='Usage: %prog [options]\n '\
                              'Clusters and catalogs FHD component arrays.')    
    parser.add_option('-r','--fhd_run',dest="fhd_run",
                      help='The fhd run name (required), e.g. \'pac_decon_eor1_June2016\'.')    
    parser.add_option('-b','--fhd_path',dest='fhd_base',
                      help='FHD output directory (optional). Defaults to \n%s'%fp.fhd_base())
    parser.add_option('-o','--obslist_file',dest='obslist_file',
                      help='File containing list of obsids to process (optional). Defaults to all obsids.')    
    parser.add_option('-e','--eps_factor',dest='eps_factor',
                      help='The clustering radius in units of beam width (optional). Defaults to 0.5.')
    parser.add_option('-m','--max_comps',dest='max_comps',
                      help='The maximum number of components to use (optional).')
    parser.add_option('-s','--shift_ra',dest='shift_ra',
                      help='The coordinate below which to shift RA +360 deg for continuity (optional). Defaults to 180.')
    parser.add_option('-f','--re_fetch_comps',action='store_true',dest='re_fetch_comps', 
                      help='Flag to re-fetch FHD component data. Default is False.')
    parser.add_option('-c','--re_cluster_comps',action='store_true',dest='re_cluster_comps',
                      help='Flag to re-cluster components. Default is False.')
    parser.add_option('-k','--re_catalog_sources',action='store_true',dest='re_catalog_sources',
                      help='Flag to re-catalog sources after clustering. Default is False.')

    (options,args) = parser.parse_args()

    return options,args
   

if __name__=="__main__":
    options,args = get_args()
   
    if options.fhd_run is None: raise ValueError,'-r --fhd_run cannot be None.'
    else: fhd_run = options.fhd_run

    if options.fhd_base is not None: fp.set_fhd_base(options.fhd_base)

    if options.obslist_file is None:  obsids = fp.get_obslist(fhd_run)
    else: obsids = np.loadtxt(options.obslist_file).astype(int)
    
    if options.eps_factor is None: eps_factor = 0.5
    else: eps_factor = float(options.eps_factor)

    if options.max_comps is not None: max_comps = int(options.max_comps)
    else: max_comps = options.max_comps

    if options.shift_ra is None: shift_ra = 180.
    else: shift_ra = float(options.shift_ra)

    re_fetch_comps = options.re_fetch_comps
    re_cluster_comps = options.re_cluster_comps
    re_catalog_sources = options.re_catalog_sources
    
    s = '%sfhd_%s'%(fp.fhd_base(),fhd_run)
    assert os.path.exists(s)
    kgs_out = '%sfhd_%s/katalogss/'%(fp.fhd_base(),fhd_run)
    if not os.path.exists(kgs_out): os.mkdir(kgs_out)
    comp_file = kgs_out+'%s_components.p'%fhd_run
    cat_file = kgs_out+'%s_catalogs.p'%fhd_run
   
    beamxx,beamyy,residual = get_images()
    obsids = remove_badobs()

    if not re_fetch_comps and os.path.exists(comp_file):
        comps,meta = pickle.load(open(comp_file))
    else: comps,meta = get_component_data()

    if not meta['clustered'] or re_cluster_comps: 
        comps,meta = cluster_components()

    if not os.path.exists(cat_file) or re_catalog_sources: 
        generate_catalog()
    
    print 'done'
