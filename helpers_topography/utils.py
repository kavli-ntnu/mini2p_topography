### GENERAL HELPERS FOR TOPOGRAPHY


####### FOR SMOOTHING FUNCTION ################################################
import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel, Gaussian2DKernel

###### FOR BINNING FUNCTIONS ##################################################
from scipy.stats import binned_statistic_2d, circmean, circstd



#### GENERAL 
def chunks(lst, n):
    ''' Generator yielding n sized successive chunks from list lst'''
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def abs_min_diff_angle(angle1, angle2):
    ''' 
    Given two angles in [-pi,pi] or [0, 2*pi] (radians!), 
    get min absolute difference in between them 

    Parameters
    ----------
    angle1 : float : radians
    angle2 : float : radians 

    Returns
    -------
    min_angle : float : minimum angle
    '''
    #a = angle2 - angle1
    #a = np.mod((a + np.pi), 2*np.pi) - np.pi

    # Make sure angle is [0, 2*pi]
    angle1 = (angle1 + 2 * np.pi) % (2 * np.pi)
    angle2 = (angle2 + 2 * np.pi) % (2 * np.pi)
    diff = np.diff([angle1, angle2])

    min_angle = np.min([2*np.pi-np.abs(diff), np.abs(diff)])
    return min_angle

def corr2(a,b):
    '''
    2-D correlation coefficient 
    Python equivalent to matlab corr2
    (check https://stackoverflow.com/questions/29481518/python-equivalent-of-matlab-corr2)
    Treats nans as zeros and fills masked numpy arrays with zeros where masked.

    '''

    if isinstance(a, np.ma.core.MaskedArray):
        a = np.ma.filled(a.copy(), fill_value=0)
    if isinstance(b, np.ma.core.MaskedArray):
        b = np.ma.filled(b.copy(), fill_value=0)

    a = np.nan_to_num(a.copy())
    b = np.nan_to_num(b.copy())
    
    def mean2(x):
        y = np.sum(x) / np.size(x);
        return y
    
    a = a - mean2(a)
    b = b - mean2(b)

    r = (a*b).sum() / np.sqrt((a*a).sum() * (b*b).sum());
    return r




#### SMOOTHING FUNCTIONS 
def smooth_stat2D(stat2D, sigma=2):
    ''' Smooth 2D histogram '''
    kernel = Gaussian2DKernel(x_stddev=sigma) # Astropy
    masked_stat2D = np.ma.masked_invalid(stat2D) 

    pad_width = int(5*sigma)
    binned_signal_padded = np.pad(masked_stat2D, pad_width=pad_width, mode='symmetric') 
    binned_signal_smoothed = convolve(binned_signal_padded, kernel, boundary='extend', preserve_nan=True)\
        [pad_width:-pad_width, pad_width:-pad_width] 
    
    binned_signal_smoothed = np.ma.masked_where(masked_stat2D.mask, binned_signal_smoothed)
    binned_signal_smoothed = binned_signal_smoothed.filled(fill_value=np.nan)
    return binned_signal_smoothed

def smooth_stat1D(stat1D, sigma=2):
    ''' Smooth 1D histogram '''
    kernel = Gaussian1DKernel(stddev=sigma) # Astropy
    masked_stat1D = np.ma.masked_invalid(stat1D) 

    pad_width = int(5*sigma)
    binned_signal_padded = np.pad(masked_stat1D, pad_width=pad_width, mode='symmetric') 
    binned_signal_smoothed = convolve(binned_signal_padded, kernel, boundary='extend', preserve_nan=True)[pad_width:-pad_width] 

    binned_signal_smoothed = np.ma.masked_where(masked_stat1D.mask, binned_signal_smoothed)
    binned_signal_smoothed = binned_signal_smoothed.filled(fill_value=np.nan)
    return binned_signal_smoothed


#### BINNINGS STATISTICS (FOV, ...) 

class BinningsStats:
    '''
    Class for binnings statistics 

    The purpose of this class is to collect commonly used statistics 
    used throughout in binned_statistic_1D and 2D

    '''
    
    def __init__(self, debug=True):
        self.debug = debug
        
    ########################################


    def nanmean(self, array):
        ''' Average array and handle all None / Nan cases ''' 
        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(np.float)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                return np.nanmean(array)
    
    def nanmedian(self, array):
        ''' Median array and handle all None / Nan cases ''' 
        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(np.float)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                return np.nanmedian(array)
            
            
    def nanmax(self, array):
        ''' Return max of array and handle all None / Nan cases ''' 
        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(np.float)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                return np.nanmax(array)
                
            
    def nanmin(self, array):
        ''' Return min of array and handle all None / Nan cases ''' 
        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(np.float)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                return np.nanmin(array)


    def nancircmean(self, array):
        ''' 
        Return angular mean (circmean) 
        of array and handle all None / Nan cases 
        ''' 
        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(np.float)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                # There should be no more nans, but anyway ... 
                return circmean(array, nan_policy='omit')


    def countnonzero(self, array):
        ''' 
        Count non-zero elements.
        Expects array with a bunch of zeros and ones and returns 
        the number of non zero elements.
        ''' 

        if not len(array):
            return np.nan
        else:
            array = array.copy().astype(int)
            array_filter = np.isfinite(array)
            if not np.sum(array_filter):
                return np.nan
            else:
                return len(array[array>0])