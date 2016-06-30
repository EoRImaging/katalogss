import numpy as np


def centroid(x, flux):
     mu = np.sum(x*flux)/np.sum(flux)
     sd = np.sqrt(np.sum(flux * (x-mu)**2)/np.sum(flux))
     return  mu,sd


def approx_stokes_i(Axx,Ayy):
     try: a = np.sqrt((Axx**2 + Ayy**2)/2.)
     except TypeError: 
          a = type(Axx)()
          a.header = Axx.header
          a.data = np.sqrt((Axx.data**2 + Ayy.data**2)/2.)
     return a


def sigma_clip(A,n_sigma,err=None, return_inds=False):
     A=np.array(A)
     if err is not None:
         w=1/err**2
         V1= np.sum(w)
         V2 = np.sum(w**2.)
         mu = np.sum(A*w)/np.sum(w)
         var =  np.sum(w*(A-mu)**2)/V1
         s2 = var / (1-V2/V1**2.) 
         sig = np.sqrt(s2)
     else: mu,sig = np.mean(A),np.std(A)
     wa=np.where(abs(A-mu)<n_sigma*sig)[0]
     if return_inds:  return [A,wa]
     else: return A[wa]


def weighted_mean(A,sig):
     w=1./sig**2
     V1= np.sum(w)
     V2 = np.sum(w**2.)
     mu = np.sum(A*w)/np.sum(w)
     sig_mu = np.sqrt(1./np.sum(w))
     var =  np.sum(w*(A-mu)**2)/V1
     s2 = var / (1-V2/V1**2.)
     sig = np.sqrt(s2)
     return [mu,sig_mu,sig]


def header_keys(file,case='None'):
    h = open(file).readline().strip('#').strip('\n').split()
    if case=='lower': h=[hi.lower() for hi in h]
    if case=='upper': h=[hi.upper() for hi in h]
    return h


def getbinsize(A):
    return 3.5*np.std(A)/len(A)**(1/3.)


def getbins(A,binsize=None):
    if binsize is None: bs=getbinsize(A)
    nbins=np.ceil((max(A)-min(A))/binsize)+1
    diff=nbins*binsize - (max(A)-min(A))
    bins=np.arange(min(A)-diff/2,max(A)+diff/2+binsize,binsize)
    return bins


def minmax(a):
    return min(a),max(a)


def span(a):
    return max(a)-min(a)

