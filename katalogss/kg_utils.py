'''
Generally useful functions.
'''

import numpy as np


def centroid(x, weights):
     '''
     Return the weighted mean and rms width.

     Parameters
     ----------
     x: `array`
         Values to be centroided.
     weights: `array`
         Values to weight `x` by. 
     '''
     mu = np.sum(x*weights)/np.sum(weights)
     sd = np.sqrt(np.sum(weights * (x-mu)**2)/np.sum(weights))
     return  mu,sd


def approx_stokes_i(Axx,Ayy):
     '''
     Return the approximate Stokes I image from two XX and YY images.

     Parameters
     ----------
     Axx: `astropy.fits.HDU` object
         XX pol image.
     Ayy: `astropy.fits.HDU` object
         YY pol image.
     '''
     try: a = np.sqrt((Axx**2 + Ayy**2)/2.)
     except TypeError: 
          a = type(Axx)()
          a.header = Axx.header
          a.data = np.sqrt((Axx.data**2 + Ayy.data**2)/2.)
     return a


def weighted_mean(A,sig):
     '''
     Apply sigma-clipping to remove outliers.

     Parameters
     ----------
     A: `array`
         Input array of values to be clipped. 
     sig: `array`
        Errors used to determine weighting (1/`err`**2).
     '''
     w=1./sig**2
     V1= np.sum(w)
     V2 = np.sum(w**2.)
     mu = np.sum(A*w)/np.sum(w)
     var =  np.sum(w*(A-mu)**2)/V1
     s2 = var / (1-V2/V1**2.)
     sig = np.sqrt(s2)
     return [mu,sig]


def sigma_clip(A,n_sigma,err=None, return_inds=False):
     '''
     Apply sigma-clipping to remove outliers.

     Parameters
     ----------
     A: `array`
         Input array of values to be clipped. 
     n_sigma: `float`
         Number of standard deviations outside which to clip.
     err: `array`, optional
        Errors used to determine weighting (1/`err`**2).
     return_inds: `bool`, optional
         If True, return `A` unmodified and indices of non-outliers. Default is False.
     '''
     A=np.array(A)
     if err is not None: mu,sig = weighted_mean(A,err)
     else: mu,sig = np.mean(A),np.std(A)
     wa=np.where(abs(A-mu)<n_sigma*sig)[0]
     if return_inds:  return [A,wa]
     else: return A[wa]



def column_names(s,delim=''):
     '''
     Return the column names from the first line of an text data file.
     
     Parameters
     ----------
     s: `file` object or `string`
        The data file from which to read the header column names (e.g. ``s = open(<filename>)``). The header line may also be passed as a string (e.g. ``s = open(<filename>).readline()``). 
     delim: `string`, optional
         Delimiter separating column names. Defaults to whitespace.
     '''
     if type(s) is file: s = s.readline()
     h = open(filename).readline().strip('#').strip('\n').split(delim)
     h=[hi.lower().strip(' ') for hi in h]
     return h


def getbinsize(A):
     '''
     Estimate appropriate binsize for histogram using Scott's rule.

     Parameters
     ----------
     A: `array`
         Values for which to estimate an appropriate bin size. 
     '''
     return 3.5*np.std(A)/len(A)**(1/3.)


def getbins(A,binsize=None):
     '''
     Return bins for a given `array` for a given (optional) bin size. If not specified, use Scott's rule. 

     Parameters
     ----------
     A: `array`
         Values for which to estimate an appropriate bin size. 
     binsize: `float`, optional
         The desired width of the bins. If None, use Scott's rule. 
     '''
     if binsize is None: bs=getbinsize(A)
     nbins=np.ceil((max(A)-min(A))/binsize)+1
     diff=nbins*binsize - (max(A)-min(A))
     bins=np.arange(min(A)-diff/2,max(A)+diff/2+binsize,binsize)
     return bins


def minmax(a):
     '''
     Return the minimum and maxium of an array.

     Parameters
     ----------
     a: `array`
         Values to return min and max of.
     '''
     return min(a),max(a)


def span(a):
     '''
     Return the difference between the max and min of an array.

     Parameters
     ----------
     a: `array`
         Values to return min and max of.
     '''
     return np.diff(minmax(a))

