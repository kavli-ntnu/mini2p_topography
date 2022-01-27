### Figure S7 - Obenhaus et al.
# panel D: anatomical masks
# Plots an overlay of anatomical masks as 
# annotated by the experimenter


import sys, os
import os.path
from pathlib import Path
from matplotlib.figure import Figure
import numpy as np 
import pandas as pd 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS #########################################################################################
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from helpers_topography.fov_stitching import make_composite_fov
from dj_schemas.anatomical_alignment import get_aligned_anatomical_masks
from skimage.exposure import rescale_intensity
from PIL import Image, ImageCms
import colorsys

rgb_to_hsv = np.vectorize(colorsys.rgb_to_hsv)
hsv_to_rgb = np.vectorize(colorsys.hsv_to_rgb)

##### LOAD SCHEMA COMPONENTS ########################################################################## 
from dj_schemas.dj_conn import * 


##### EXPORT LOCATION #################################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'
figure_dir = Path(figure_dir)


#### HELPERS FOR CMYK -> RGB conversion ###############################################################
def _convert_to_8bit(data):
    '''
    Convert numpy array from arbitrary range
    to np.unit8 
    '''
    data = np.nan_to_num(data)
    data = rescale_intensity(data, in_range='image', out_range=np.uint8)
    data = data.astype(np.uint8)
    return  np.nan_to_num(data)

def _shift_hue(arr, 
               hout, 
               ):
    r, g, b, a = np.rollaxis(arr, axis=-1)
    h, s, v = rgb_to_hsv(r, g, b)
    h += hout
    r, g, b = hsv_to_rgb(h, s, v)
    arr = np.dstack((r, g, b, a))
    return arr

def _colorize(image, 
              hue_change, 
              ):
    '''
    Hue change PIL image `image` with the given
    `hue_change` (hue within 0-360)
    returns another PIL image

    '''
    img = image.convert('RGBA')
    arr = np.array(np.asarray(img).astype('float'))
    hue_shifted = _shift_hue(arr, 
                             hue_change/360., 
                             )
    new_img = Image.fromarray(hue_shifted.astype('uint8'), 'RGBA')

    return new_img

#######################################################################################################

def _plot_fov_and_masks(animal_name, 
                        projection_short, 
                        ssim_thresh, 
                        center_plane=0,
                        plot_mask_fov='fov',
                        percentile=99.8,
                        size_scalebar=50.,
                        ):
    '''
    Anatomical FOV + 
    PLOTTING OF ANATOMICAL MASKS! 


    '''
    assert plot_mask_fov in ['fov','mask'], f'Choose "fov" or "mask" for plot_mask_fov'

    # Retrieve sessions and plot FOV composite 
    session_filter = (Session.proj(..., metasession_ref='metasession_name')
                                   & f'animal_name = "{animal_name}"'
                                   & FilteredSessions)
                      


    metasessions   = (AlignmentFOV.proj(diff_ssim = 'ssim_warped-ssim_original') 
                                  & session_filter 
                                  & f'diff_ssim > {ssim_thresh}'
                                  & f'projection_short = "{projection_short}"')
                     

    alignment_composite, x_min, x_max, y_max, y_min  = make_composite_fov(metasessions, 
                                                                          center_plane = center_plane, 
                                                                          padding_composite = 5.,
                                                                          add_edge=False,
                                                                          )

    if plot_mask_fov == 'mask':

        # Retrieve anatomical masks (manual annotations)
        region_masks_mec = get_aligned_anatomical_masks(animal_name=animal_name,
                                                region='mec_label',
                                                projection_short=projection_short,
                                                edge=False
                                                )
        region_masks_pas = get_aligned_anatomical_masks(animal_name=animal_name,
                                                    region='pas_label',
                                                    projection_short=projection_short,
                                                    edge=False
                                                    )
        # Normalize
        region_masks_mec /= region_masks_mec.max()
        region_masks_pas /=region_masks_pas.max()

        # Create a CMYK stack from those masks, 
        # filling "C", and "Y" and leaving the others at zero
        # See also my gist: 
        # https://gist.github.com/horsto/072758ba4b6a292264cd65a9a98a80d2

        all_masks_cmyk = np.stack([region_masks_mec,  
                                   np.zeros_like(region_masks_mec), 
                                   region_masks_pas, 
                                   np.zeros_like(region_masks_mec)
                                  ],
                                  axis=2)
        all_masks_cmyk8bit = _convert_to_8bit(all_masks_cmyk)

        # profiles
        cmyk_profile = '/Users/hotte/Documents/python/hotte-dj-moser-imaging/color_profiles/USWebCoatedSWOP.icc'
        rgb_profile  = '/Users/hotte/Documents/python/hotte-dj-moser-imaging/color_profiles/sRGB Color Space Profile.icm'

        img_pil = Image.fromarray(all_masks_cmyk8bit, mode='CMYK')
        img = ImageCms.profileToProfile(img_pil, 
                                      cmyk_profile, 
                                      rgb_profile, 
                                      renderingIntent=0, 
                                      outputMode='RGB')
        # HUE CHANGE
        img_hue = _colorize(img, hue_change=185)   
        # ... convert back to RGB image      
        all_masks_rgb = np.array(img_hue)




    # Plot ...  
    cmap = 'Greys'
    figure = plt.figure(figsize=(10,10))
    ax = figure.add_subplot(111)

    if plot_mask_fov == 'fov':
        ax.imshow(alignment_composite, vmin=0, vmax=np.nanpercentile(alignment_composite, percentile), cmap=cmap)
    else: 
        ax.imshow(all_masks_rgb)

    # Scalebar ... 
    if size_scalebar is not None:
        ax.plot([x_max-size_scalebar-10,x_max-size_scalebar-10+size_scalebar], 
                [y_max-10,y_max-10], 
                lw=3, 
                color='k', 
                solid_capstyle='butt'
                )

    ax.set_xlim([x_min,x_max]), ax.set_ylim([y_max,y_min])

    # Construct title 
    projection_key = (metasessions * FOVProjParam).fetch('projection_key')[0]
    ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([])    
    sns.despine(left=True, bottom=True)

    print(f'Animal {animal_name} (n={len(metasessions)+1} sessions) {projection_key}')
    print('Saving figure under {}'.format(str(figure_dir)))
    figure.savefig(figure_dir / f'{animal_name} composite {plot_mask_fov}.png', bbox_inches='tight', dpi=300)



if __name__ == "__main__":

    #### PLOT COMPOSITE FOV EXAMPLE #####################################################
    ssim_thresh = 0. # This is: [ssim_warped-ssim_original] in AlignmentFOV
    center_plane = 0
    animal_name = '88592'
    projection_short = 'mean_im2' # FOVProjParam()
    percentile = 99.8

    mask_fov = 'fov' # 'mask' or 'fov' - either plot the anatomical fov defined by projection_short
                     # or the annotated region masks for MEC / PAS only
    size_scalebar = 50. 


    print('Plotting (anatomical) FOV composite')
    _plot_fov_and_masks(animal_name, 
                        projection_short, 
                        ssim_thresh, 
                        center_plane=0,
                        plot_mask_fov=mask_fov,
                        percentile=percentile,
                        size_scalebar=size_scalebar,
                        )
   

    print('Success.')