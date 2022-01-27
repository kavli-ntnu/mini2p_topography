### HELPER FUNCTIONS FOR FOV STITCHING / EVALUATION
from pathlib import Path
import numpy as np 
import pandas as pd 
from statannot import add_stat_annotation
from skimage.transform import EuclideanTransform, AffineTransform, warp
from skimage.metrics import structural_similarity, mean_squared_error
from skimage.feature import canny
from scipy.stats import binned_statistic_2d


#### PLOTTING ################################################################################
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Load more colormaps 
import cmasher as cmr

####### LOAD SCHEMATA ########################################################################
import datajoint as dj
from dj_schemas.dj_conn import * 

# UTILS 
from .utils import BinningsStats, smooth_stat2D
from dj_schemas.constants import COMPRESSION_CUTOFF, SNR_CUTOFF


# Helper for constructing data dictionaries 


def get_cell_data(cells, table_params, additional_columns=None):
    ''' 
    Retrieve x, y and score info for cells (DJ join object)  
    
    Attributes
    ----------
    cells : datajoint join object containing all necessary columns 
    table_params : dictionary 
                   'Score', 
                   'OtherScore', (can be None)
                   'attributes'  (can be None)
    additional_columns : list
                    Additional column names to add to 
                    pandas dataframe

    Returns 
    -------
    cells_data : pandas dataframe containing fetched cell data


    Note: 
    This is used both in the process of score neighbourhood analysis,
    as well as the macroscopic (FOV) based analysis 

    '''
    assert isinstance(table_params, dict)
    Score      = table_params['Score']
    OtherScore = table_params['OtherScore']
    attributes = table_params['attributes']

    assert (Score in cells.heading.names), f'{Score} not in provided join object'
   
    final_columns_df = ['cell_id', 'center_x_corr', 'center_y_corr', Score]

    if additional_columns is not None:
        assert isinstance(additional_columns, list)
        final_columns_df.extend(additional_columns)

    # Check if "OtherScore" is there 
    if OtherScore is not None: 
        if not hasattr(OtherScore):
            # ... i.e. that is a real score and not a calculated one
            assert (OtherScore in cells.heading.names), f'{OtherScore} not in provided join object' 
            final_columns_df.append(OtherScore)

    # Check if "attributes" are there 
    if attributes is not None: 
        for attr in attributes:
            assert (attr in cells.heading.names), f'{attr} not found in join object'
            final_columns_df.append(attr)

    final_columns_df = list(set(final_columns_df))

    cells_data = pd.DataFrame(cells.fetch(as_dict=True))
    cells_data = cells_data[final_columns_df]
    cells_data.set_index('cell_id', inplace=True)
    return cells_data



# Helpers for expanding FOV 
def shift_image(image, vector):
    ''' Shift image by vector. Vector is a tuple of x / y shifts '''
    assert type(vector) == tuple, 'Vector is not tuple'
    assert len(vector) == 2, 'Vector has to be tuple of length 2'
    
    
    transform = EuclideanTransform(translation=vector)
    shifted = warp(image, transform, mode='constant', preserve_range=True)
    shifted = shifted.astype(image.dtype)
    return shifted

def shift_points(points, vector):
    ''' Shift points by vector. Vector is a tuple of x / y shifts '''
    assert type(vector) == tuple, 'Vector is not tuple'
    assert len(vector) == 2, 'Vector has to be tuple of length 2'
    
    points = points.copy()
    shifted = []
    for point in points:
        point[0] -= vector[0]
        point[1] -= vector[1]
        shifted.append(point)
    return np.stack(shifted)


def align_session_imgs(projection, 
                       points_left, 
                       points_right, 
                       trans_ops, 
                       transform=False, 
                       transform_method='affine'
                       ):
    '''
    Helper for shifting projection and alignment points 
    and building EuclideanTransform / AffineTransform / ... object

    if "transform == True":
        Build transform object and return scikit image transformer and 
        warped image

    '''
    # First pad
    projection_padded  = (np.pad(projection, trans_ops['padding']))
    # ... then shift 
    projection_shifted  = shift_image(projection_padded, trans_ops['vector'])

    if transform: 
        # Shift points 
        points_left_shifted      = shift_points(points_left,  trans_ops['vector_corr'])
        points_right_shifted     = shift_points(points_right, trans_ops['vector_corr'])
        
        # Choose transform method
        assert transform_method in ['affine','euclidean'], f'Method "{transform_method}" not known'
        if transform_method == 'affine':
            tform = AffineTransform()
        elif transform_method == 'euclidean':
            tform = EuclideanTransform()

        status = tform.estimate(points_left_shifted, points_right_shifted)
        if not status:
            raise ValueError('{transform_method} transform ran into an error during fitting')      
            
        # Unwarp image 
        picture_unwarped = warp(projection_shifted, tform, order=4, mode='constant',\
                                preserve_range=False, clip=False)
    else:
        tform = None
        picture_unwarped = None
        
    return projection_shifted, picture_unwarped, tform



# Helper for extracting boundaries of image 
def get_x_y_lim(composite, thresh=1 , padding=20):
    composite = np.nan_to_num(composite)
    meanx = np.mean(composite,axis=0)
    crossingx_1 = np.where(meanx>thresh)[0][0] # x lim min
    crossingx_2 = np.where(meanx>thresh)[0][-1] # x lim max
    meany = np.mean(composite,axis=1)
    crossingy_2 = np.where(meany>thresh)[0][0] # y lim min
    crossingy_1 = np.where(meany>thresh)[0][-1] # y lim max
    return crossingx_1-padding, crossingx_2+padding, crossingy_1+padding, crossingy_2-padding

def crop_center(img,cropx,cropy):
    if (img.shape[1] > cropx) or (img.shape[0] > cropy): 
        print('Image to be center cropped is not big enough. Padding.')
    img = np.pad(img.copy(), 500, mode='constant', constant_values=0) # Pad with arbitrarily large value

    # Center crop
    y,x = img.shape
    startx = x//2-(cropx//2)
    starty = y//2-(cropy//2)    
    return img[starty:starty+cropy,startx:startx+cropx]


# Calculate overlapping coordinates
def get_center_crop(x_min_unwarped, 
                    x_max_unwarped, 
                    y_min_unwarped, 
                    y_max_unwarped, 
                    x_min_reference, 
                    x_max_reference, 
                    y_min_reference, 
                    y_max_reference):

    ''' 
    Get max of the two minimal x coordinates, 
    min of the two maximal x coordinates, etc... 
    to extract image core (where they overlap maximally)
    '''
    return np.max([x_min_unwarped, x_min_reference]), \
           np.min([x_max_unwarped, x_max_reference]), \
           np.min([y_min_unwarped, y_min_reference]), \
           np.max([y_max_unwarped, y_max_reference])


def plot_metrics(metrics_df, 
                 plot_ssim=True, 
                 plot_mse=False, 
                 show_stats=True, 
                 ssim_thresh=None, 
                 mse_thresh=None, 
                 marker_size=9,
                 ax_ssim=None,
                 ax_mse=None,
                 save_path=None):

    ''' Plot SSIM and MSE metrics from comparisons of warped and original pictures '''
    
    # Set up a grid to plot survival probability against several variables
    sns.set(font_scale=1.5, style='white')
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure_for_export = []
    ###### SSIM ########################################################################################################################
    if plot_ssim:
        if ax_ssim is None:
            figure_ssim = plt.figure(figsize=(1.5,5))
        
        sns.stripplot(data=metrics_df[['ssim_original','ssim_warped']],
                      s=marker_size,
                      linewidth=0, 
                      alpha=.75, 
                      zorder=0, 
                      palette='Greys',
                      ax=ax_ssim)
        g = sns.pointplot(data=metrics_df[['ssim_original','ssim_warped']], 
                          ci='sd', 
                          color='k', 
                          lw=1, 
                          zorder=10,
                          ax=ax_ssim)

        sns.despine()
        g.set_xticklabels(['Original', 'Aligned'])
        g.set_xticklabels(g.get_xticklabels(), rotation=45, ha='right', va='top')
        
        #yticks = g.get_yticks()    
        yticks_ssim = [0, .2, .4, .6, .8]
        g.set_yticks(yticks_ssim)
        g.set_yticklabels(yticks_ssim)

        plt.gca().margins(x=0.2)

        g.set_ylabel('SSIM')
        if show_stats:
            add_stat_annotation(g, data=metrics_df, box_pairs=[('ssim_original','ssim_warped')], 
            test='Mann-Whitney', verbose=0, loc='inside')
        if ssim_thresh is not None:
            g.axhline(y=ssim_thresh, ls=':', color='k')

        sns.despine(left=True, bottom=True)
        if save_path is not None: 
            save_path = Path(save_path)
            print('Saving figure under {}'.format(str(save_path)))
            figure_ssim.savefig(save_path / f'SSIM alignment metrics (n={len(metrics_df)}).pdf', bbox_inches='tight')
        

        if ax_ssim is None:
            figure_for_export.append(figure_ssim)
        
    ###### MSE  ########################################################################################################################
    if plot_mse:
        if ax_mse is None:
            figure_mse = plt.figure(figsize=(1.5,5))
        
        sns.stripplot(data=metrics_df[['mse_original','mse_warped']],
                      s=marker_size,
                      linewidth=0, 
                      alpha=.75, 
                      zorder=0, 
                      palette='Greys',
                      ax=ax_mse)
        g = sns.pointplot(data=metrics_df[['mse_original','mse_warped']], 
                          ci='sd', 
                          color='k', 
                          lw=1, 
                          zorder=10,
                          ax=ax_mse)

        sns.despine(left=True, bottom=True)
        g.set_xticklabels(['Original', 'Aligned'])
        g.set_xticklabels(g.get_xticklabels(), rotation=45,  ha='right', va='top')

        yticks_mse = g.get_yticks()    
        g.set_yticklabels(yticks_mse)

        plt.gca().margins(x=0.2)

        g.set_ylabel('MSE')
        if show_stats:
            add_stat_annotation(g,data=metrics_df, box_pairs=[('mse_original','mse_warped')], 
            test='t-test_paired', verbose=0,loc='inside')
        if mse_thresh is not None:
            g.axhline(y=mse_thresh, ls=':', color='k')

        if save_path is not None: 
            save_path = Path(save_path)
            print('Saving figure under {}'.format(str(save_path)))
            figure_mse.savefig(save_path / f'MSE alignment metrics (n={len(metrics_df)}).pdf', bbox_inches='tight')
        

        if ax_mse is None:
            figure_for_export.append(figure_mse)

    return figure_for_export



def evaluate_unwarping(picture_unwarped, picture_reference, projection_left, projection_right, \
                                y_range_left, x_range_left, y_range_right, x_range_right,  verbose=True):
    '''
    Evaluate success of alignment/unwarping by comparing pictures (picture_reference, picture_unwarped) to their 
    non-warped (raw) counterparts (projection_left, projection_right).
    Crop a center region to compare only "valid" areas.

    Note: The naming "left"/"right" stems from the GUI used in the process of defining anatomical landmarks.

    Returns SSIM and mean squared error between those core regions. 
    
    Parameters
    ----------
    picture_unwarped : np.array (2D)
    picture_reference: np.array (2D)
    projection_left:   np.array (2D)
    projection_right:  np.array (2D)

    ranges (y_... x_...): float

    Returns
    -------
    SSIM original : structural similarity index (SSIM) of raw (non-warped) pictures
    SSIM warped   :                ...                 of warped pictures
    MSE original  : mean squared error between raw (non-warped) pictures
    MSE warped    : mean squared error between warped pictures
    '''


    try:
        x_min_unwarped, x_max_unwarped, y_min_unwarped, y_max_unwarped     = get_x_y_lim(picture_unwarped, padding=0)
        x_min_reference, x_max_reference, y_min_reference, y_max_reference = get_x_y_lim(picture_reference,padding=0)
    except IndexError:
        print('Boundaries of image could not be extracted (possbibly because it is empty)')
        return np.nan, np.nan, np.nan, np.nan 

    # Crop valid center    
    x_min_overlap, x_max_overlap, y_min_overlap, y_max_overlap = get_center_crop(
                                                                            x_min_unwarped, 
                                                                            x_max_unwarped, 
                                                                            y_min_unwarped, 
                                                                            y_max_unwarped, 
                                                                            x_min_reference, 
                                                                            x_max_reference, 
                                                                            y_min_reference, 
                                                                            y_max_reference)
    
    # ... crop to core region
    unwarped_core  = picture_unwarped[y_max_overlap:y_min_overlap, x_min_overlap:x_max_overlap]
    reference_core = picture_reference[y_max_overlap:y_min_overlap, x_min_overlap:x_max_overlap]
    
    # ... also crop original pictures to core
    projection_left_cropped  = projection_left[y_range_left[0]:y_range_left[1],x_range_left[0]:x_range_left[1]]
    projection_right_cropped = projection_right[y_range_right[0]:y_range_right[1],x_range_right[0]:x_range_right[1]]
    cropped_left  = crop_center(projection_left_cropped,unwarped_core.shape[1],unwarped_core.shape[0])
    cropped_right = crop_center(projection_right_cropped,unwarped_core.shape[1],unwarped_core.shape[0])
    
    # Original and warped core pictures should have the same dimensions
    assert cropped_left.shape == cropped_right.shape == unwarped_core.shape == reference_core.shape, 'Pictures have different dimensions'
    
    # Calculate similarity metrics
    ssim_original = structural_similarity(cropped_left, cropped_right)
    ssim_warped = structural_similarity(unwarped_core, reference_core)
    
    mse_original = mean_squared_error(cropped_left, cropped_right)
    mse_warped   = mean_squared_error(unwarped_core, reference_core)
    
    if verbose: print('[warped/original] SSIM diff: {:.2f} | MSE ratio: {:.2f}'.format(ssim_warped-ssim_original,mse_warped/mse_original))
        
    return ssim_original, ssim_warped, mse_original, mse_warped


def make_composite_fov(metasessions, 
                       center_plane=-1, 
                       padding_composite=50,
                       add_edge=True,
                       edge_value=6000.):
    '''
    Generate stitched FOV composite 
    
    metasessions : datajoint object with metasessions to be looked up in AlignmentFOV
    center_plane : plane number, defaults to -1: Last one
    padding_composite : get_x_y_lim input
    add_edge : boolean : Add a thing edge around the outline of each projection image? 
    edge_value : float : If 'add_edge' == True: Which value (grey scale) should that border have?
    '''

    # Fix threshold for edge detection(s) for now
    thresh = 1

    def __add_edge(image, thresh=1):
        image = image.copy()

        image_mask = image.copy()
        image_mask[image_mask==0] = np.min(image_mask)
        image_mask[image_mask>thresh] = image.max()
        image_mask[image_mask<thresh] = 0

        edges = canny(
            image=image_mask,
            sigma=4,
            low_threshold=0.,
            high_threshold=.99999,
            mask=None,
            use_quantiles=True,
        )
        image += (edges.astype(int) * edge_value)
        return image


    alignment_composite = None
    for metasess in metasessions.proj(): 
        if alignment_composite is None:
            image = (AlignmentFOV & metasess).fetch1('padded_ref') # [center_plane].copy()
            if isinstance(image, dict):
                if center_plane == -1: # Find out the last one 
                    center_plane_ = list(image.keys())[center_plane]
                else:
                    center_plane_ = center_plane
                image = image[center_plane_]
                if add_edge: image = __add_edge(image, thresh=thresh)
            alignment_composite = image.copy()

        # ... if not the first result, retrieve the "aligned" image
        image = (AlignmentFOV & metasess).fetch1('padded_align')
        if isinstance(image, dict):
            if center_plane == -1: # Find out the last one 
                center_plane_ = list(image.keys())[center_plane]
            else:
                center_plane_ = center_plane
            image = image[center_plane_]
            if add_edge: image = __add_edge(image, thresh=thresh)
        alignment_composite += image.copy()

    # ... and subsequently extract x / y limits
    x_min, x_max, y_max, y_min = get_x_y_lim(alignment_composite, padding=padding_composite, thresh=thresh)
    return alignment_composite, x_min, x_max, y_max, y_min 





def plot_composite_fov(metasessions, 
                       animal_name, 
                       center_plane=0,
                       padding_composite=30, 
                       cmap='Greys', 
                       percentile=98.0, 
                       size_scalebar=50,
                       return_bounds=True, 
                       add_edge=True,
                       edge_value=4000,
                       save_path=None
                       ):
    '''
    Convenience function that combines "make_composite_fov" (just above) and 
    and creates / displays figure 
    
    '''
    from dj_schemas.anatomical_alignment import FOVProjParam
    alignment_composite, x_min, x_max, y_max, y_min  = make_composite_fov(metasessions, 
                                                                          center_plane = center_plane, 
                                                                          padding_composite = padding_composite,
                                                                          add_edge=add_edge,
                                                                          edge_value=edge_value,
                                                                          )

    # Show ...  
    figure = plt.figure(figsize=(10,10))
    plt.imshow(alignment_composite, vmin=0, vmax=np.nanpercentile(alignment_composite, percentile), cmap=cmap)
    plt.plot([x_max-size_scalebar-10,x_max-size_scalebar-10+size_scalebar], [y_max-10,y_max-10], lw=3, color='k', solid_capstyle='butt')

    plt.xlim([x_min,x_max]), plt.ylim([y_max,y_min])


    # Construct title 
    projection_key = (metasessions * FOVProjParam).fetch('projection_key')[0]
    plt.title(f'Animal {animal_name} (n={len(metasessions)+1} sessions) {projection_key}', y=.99)
    plt.gca().get_xaxis().set_ticks([]);plt.gca().get_yaxis().set_ticks([])    
    sns.despine(left=True, bottom=True)




    if save_path is not None: 
        save_path = Path(save_path)
        print('Saving figure under {}'.format(str(save_path)))
        figure.savefig(save_path, bbox_inches='tight', dpi=300)
    #plt.show()

    if return_bounds: 
        return x_min, x_max, y_max, y_min


######################################################################################################################################################
######################################################################################################################################################

##### ROIS 

def shift_px_coord(df, vector_corr):
    ''' Shift x/y pixel positions of cells (ROIs)'''
    
    x_y_nonshifted_all = []
    x_y_shifted_all = []
    
    for row in df.iterrows():
        x = row[1]['xpix_corr']
        y = row[1]['ypix_corr']
        # non-shifted
        x_y = np.stack([x,y]).T
        x_y_nonshifted_all.append(x_y)
        # shifted
        x_y_shifted = shift_points(x_y, vector_corr)
        x_y_shifted_all.append(x_y_shifted)
        
    return x_y_shifted_all, x_y_nonshifted_all


def unwrap_roi_x_y_value(cells_x_y, mappable_values):
    ''' Create arrays of flattened x and y coordinates, and corresponding scores, for every cell   '''
    assert len(cells_x_y) == len(mappable_values), f'Length of cells ({len(cells_x_y)}) and x/y coordinates ({len(mappable_values)}) do not match'
    mappable_values_all = []
    xs_all              = []
    ys_all              = []
    # Loop over cells and collect coords and values
    for no in range(len(cells_x_y)):
        for point in cells_x_y[no]:
            xs_all.append(point[1])
            ys_all.append(point[0])
            mappable_values_all.append(mappable_values[no])
    return np.array(xs_all), np.array(ys_all), np.array(mappable_values_all)



######################################################################################################################################################
######################################################################################################################################################

### SCORE MAP OPERATIONS ########


def get_fov_score_map(cells, params, vector_corr, tform=None, shuffle_iter=0):
    '''
    1. Get cell data (xpix_corr, ypix_corr values and score) 
    2. Shift px coordinates like previously done for all other vectors (Align ...)
    3. Unwrap and
    4. bin 2D 

    Creates shuffled distribution by scrambling relationship cell id vs. scores 


    Parameters
    ----------
    cells : datajoint join object containing all necessary info (x, y, score)
    params : dict
             - Score
             - statistic for binning 
                -> implemented as function in helpers_topography.utils BinningsStats()
             - bins_x, bins_y (binning 2D)
                   

    '''
    # Make sure we only fetched from one unique session and that there are no duplicate cell IDs
    assert len(set(cells.fetch('session_name'))) == 1
    assert len(cells.fetch('cell_id')) == len(set(cells.fetch('cell_id')))
    
    # Get cell data dataframe (score, x, y)
    cell_data = get_cell_data(cells, params, additional_columns=['xpix_corr','ypix_corr'])
    score_values = cell_data[params['Score']].values

    # Shift x/y coordinates
    x_y_shifted, _  = shift_px_coord(cell_data, vector_corr)
    
    # Check if (affine/euclidean) transform was loaded 
    if tform is not None:
        assert isinstance(tform, (EuclideanTransform, AffineTransform)) , \
            f'Transform object (tform) is of unknown type ({type(tform)})'
        print(f'Using transform on {len(x_y_shifted)} points')
        warped_x_y = [tform.inverse(x_y) for x_y in x_y_shifted]
        x_y_shifted = warped_x_y # assign back to original
        
    ys_unwrapped, xs_unwrapped, score_unwrapped = unwrap_roi_x_y_value(x_y_shifted, score_values)
    # Bin + smooth
    stat2D, y_edges, x_edges , _ = binned_statistic_2d(
                                                    x = ys_unwrapped,
                                                    y = xs_unwrapped,\
                                                    values = score_unwrapped, \
                                                    statistic = params['statistic'], 
                                                    bins = [params['bins_y'], params['bins_x']]
                                                      )

    stat2D_shuff = []
    if shuffle_iter: 
        score_values_shuff = score_values.copy()
        for _ in range(shuffle_iter): # np.min([len(score_values), shuffle_iter])
            np.random.shuffle(score_values_shuff) # This is in-place!
            ys_unwrapped_, xs_unwrapped_, score_unwrapped_ = unwrap_roi_x_y_value(x_y_shifted, score_values_shuff)
    
            # Bin + smooth
            stat2D_ , _, _, _ = binned_statistic_2d(
                                                x = ys_unwrapped_,
                                                y = xs_unwrapped_,\
                                                values = score_unwrapped_, \
                                                statistic = params['statistic'], 
                                                bins = [params['bins_y'], params['bins_x']]
                                                )
            stat2D_shuff.append(stat2D_)



    return stat2D, stat2D_shuff, x_edges, y_edges


def stack_project_fov_score_maps(scoremap_dict, statistic='nanmean', sigma=0.):
    '''
    Loop over scoremap dictionary and stack + project 
    Statistic for projection is taken from topography.utils BinningStats.
    
    '''
    ### INITIALIZE STATS #### 
    binnings_stats = BinningsStats()

    # Stack
    stacked_stat2D = []
    for scoremap in scoremap_dict.values():
        stacked_stat2D.append(scoremap)
    stacked_stat2D = np.stack(stacked_stat2D)
    
    # Project
    projection = np.apply_along_axis(getattr(binnings_stats, statistic), axis=0, arr=stacked_stat2D)
    
    if sigma > 0.:
        projection = np.nan_to_num(projection)
        projection = smooth_stat2D(projection, sigma)
    
    # Get xlim ylim 
    x_min, x_max, y_min, y_max = get_x_y_lim(np.nan_to_num(projection), thresh=-0.1, padding=0)
    
    return projection, x_min, x_max, y_min, y_max


###### BINNING PARAMETERS ########################################################################################################


class BinningParamsFOV:
    '''
    Class for binnings parameters for stitched FOV score maps 

    The purpose of this class is to fix the binnings parameters 
    used in "FOV Score maps" notebook. 

    The basic parameter setup follows the one in the main schema 
    for "get_cell_data" (helpers_topography.distributions.get_cell_data) 


    '''
    
    def __init__(self, bins_x, bins_y, debug=True, **kwargs):
        self.debug = debug

        # Load stats module from helpers_topography/distributions/get_binned_scores1D
        self.stats = BinningsStats()

        self.bins_x = bins_x
        self.bins_y = bins_y 

        self.cmap  = kwargs.get('cmap', None)

    ########################################

    @property 
    def border_bvs(self): 
        ''' 
        # Cutoffs bvs / border score
        # bvs  95th bv_field_dect_method = 'bvs' : .56
        # bvs  99th bv_field_dect_method = 'bvs' : .66
        # bvs  95th bv_field_dect_method = 'opexebo' : .68
        # bvs  99th bv_field_dect_method = 'opexebo' : .75
        # borderscore 95th : .37
        # borderscore 99th : .48

        '''
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'bvs',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (BorderScore * BVScore * CutoffsBorderScore & 'bvfield_params_id="A"' & 'bv_field_dect_method = "bvs"').proj(
                                        filtered_value='if(borderscore > borderscore_95, bvs, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Reds' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y,
                'vmin_vmax'    : (0,.7),
                }
        return params

    @property
    def grid(self):
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'gridscore',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (GridScore * CutoffsGridscore).proj(
                                        filtered_value='if(gridscore > gridscore_95, gridscore, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Reds' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y,
                'vmin_vmax'    : (0,1.3),
                }

        return params 

    @property
    def grid_99(self):
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'gridscore',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (GridScore * CutoffsGridscore).proj(
                                        filtered_value='if(gridscore > gridscore_99, gridscore, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Reds' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y
                }

        return params 


    @property
    def mvl(self):
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'mvl',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (AngularRate.Stats * CutoffsMVL).proj(
                                        filtered_value='if(mvl > mvl_95, mvl, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Blues' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y,
                'vmin_vmax'    : (0,.7),
                }

        return params 


    @property
    def mvl_cellbycell(self):
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'mvl',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (AngularRate.Stats * Shuffled.AngularRateStats).proj(
                                        filtered_value='if(mvl > mvl_95, mvl, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Blues' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y
                }

        return params 


    @property
    def ovc(self):
        params = { 
            'Score'        : 'filtered_value',
            'score_key'    : 'ovc',
            'OtherScore'   : None,
            'attributes'   : None,
            'scoretables'  : OVC.proj(
                                session_name='base_session',
                                tracking_= 'tracking_dataset', 
                                signal_= 'signal_dataset', 
                                is_ovc = 'is_ovc',
                                ovscore = 'ovscore',
                                filtered_value='if(is_ovc > 0, ovscore, 0)'
                                ),
            'score_cutoff' : None,    
            'statistic'    : self.stats.nanmax,
            'colormap'     : ['Greens' if self.cmap is None else self.cmap][0],
            'bins_x'       : self.bins_x,
            'bins_y'       : self.bins_y,
            'vmin_vmax'    : (0,.75),
            }

        return params 

    @property
    def info_content(self):
        params = {
                'Score'        : 'filtered_value',
                'score_key'    : 'information_content',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : (Ratemap.Stats * CutoffsInfoContent).proj(
                                        filtered_value='if(information_content > info_content_95, information_content, 0)'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmax,
                'colormap'     : ['Blues' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y,
                'vmin_vmax'    : (0,2.75),
                }

        return params 


    @property
    def cell_size_s2p(self):
        params = {
                'Score'        : 'npix_corr',
                'score_key'    : 'npix_corr',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : RoisCorr.proj('npix_corr','compress_ratio') * SNR.proj(snr='snr_df_f') \
                                            & 'snr > {}'.format(SNR_CUTOFF) \
                                            & 'compress_ratio > {}'.format(COMPRESSION_CUTOFF),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmean,
                'colormap'     : ['cmr.guppy_r' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y
                }

        return params 



    @property
    def snr(self):
        params = {
                'Score'        : 'snr_df_f',
                'score_key'    : 'snr_df_f',
                'OtherScore'   : None,
                'attributes'   : None,
                'scoretables'  : SNR.proj('snr_df_f'),
                'score_cutoff' : None,    
                'statistic'    : self.stats.nanmean,
                'colormap'     : ['cmr.sunburst' if self.cmap is None else self.cmap][0],
                'bins_x'       : self.bins_x,
                'bins_y'       : self.bins_y
                }

        return params 


