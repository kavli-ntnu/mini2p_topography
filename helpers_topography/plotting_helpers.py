### PLOTTING HELPERS
import math
from pathlib import Path
import numpy as np
from tqdm import tqdm, tqdm_notebook
from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# For colorbar axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import seaborn as sns

from skimage.transform import rescale
from scipy.ndimage.filters import gaussian_filter

sns.set(style='white', font_scale=1.1)

from helpers_topography.fov_stitching import (BinningParamsFOV, 
                                              stack_project_fov_score_maps)
from dj_schemas.anatomical_alignment import (ScoremapFOV, 
                                             ScoremapFOVProjParams,
                                             get_scoremaps_moran, 
                                             get_scoremaps_a_b, 
                                             prune_keys, 
                                             get_aligned_anatomical_masks)

#############################################################################################################################
#############################################################################################################################
### ANATOMICAL SCORE DISTRIBUTIONS ETC .
######### PLOTTING ########################

def crop_center(img,cropx,cropy):
    y,x = img.shape
    startx = x//2-(cropx//2)
    starty = y//2-(cropy//2)    
    return img[starty:starty+cropy,startx:startx+cropx]


#######################################################################################################################################################
#######################################################################################################################################################
# FOV maps / stats


def plot_summary_scoremaps(animal_name, 
                           binning_parameters=['grid','mvl','ovc','grid_module_1'], 
                           proj_param_id='B',
                           metasession_ref=None,
                           param_hash_id_cell='standard',
                           param_hash_session='cf83e1357eefb8bd', 
                           percentile=100.,
                           save_path=None,
                           **kwargs
                           ):
    '''
    Create summary plot of scoremaps and matching Moran's I results
    Uses:
    dj_schemas.anatomical_alignment.get_scoremaps_moran() 
    to retrieve scoremap and Moran's I entries

    Parameters
    ----------
    animal_name       :  str
                         animal_name
    binning_parameters:  list
                         List of binning parameters
                         Secondary key(binning_param_set) of ScoremapFOVParams()
    proj_param_id     :  str
                         Primary key of ScoremapFOVProjParams() 
                         Determinees projection statistic (e.g. 'nanmean') 
                         and smoothing sigma
    metasession_ref   :  str
                         Reference metasession (metasession_name).
                         For filtering if there is more than one reference session 
                         for this animal
    param_hash_id_cell : str
                         Primary key of FilteredCellsParams()
                         Default: 'standard'
    param_hash_session : str
                         Primary key of FilteredSessionsParams()   
                         Default: 'cf83e1357eefb8bd'
    percentile        :  float
                         Scaling (top percentile) for scoremap plot
    save_path         :  str or pathlib.Path
                         Save figure under this path

    **kwargs 
        bin_size_microns : float
                           Bin size in microns. If given, scale bar will be displayed.
                           Defaults to None.
        scalebar_length  : float
                           Scale bar length in microns. CAVE will only be drawn
                           if bin_size_microns is not None.
                           Defaults to 100.
        anat_mask_label  : string
                           region label as in AnatomicalMaskParams(),
                           e.g. 'mec_label' or 'pas_label'
                           If not None (default), will draw a shaded, smoothed 
                           anatomical mask on top of scoremap projection 
        anat_mask_sigma  : float
                           Smoothing sigma (Gaussian kernel) for anatomical mask 
                           Default: 1.5

    '''
    bin_size_microns    = kwargs.get('bin_size_microns', None)
    scalebar_length     = kwargs.get('scalebar_length', 100.)
    anat_mask_label     = kwargs.get('anat_mask_label', None)
    anat_mask_sigma     = kwargs.get('anat_mask_sigma', 1.5)


    if save_path is not None: 
        if isinstance(save_path, str):
            save_path = Path(save_path) 

    # Retrieve binning parameter set 
    # This is used exclusively to obtain colormaps and vmin vmax entries for display
    binning_params = BinningParamsFOV(None, None)

    sns.set(style='white', font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    figure = plt.figure(constrained_layout=False, figsize=(45,10))
    gs = figure.add_gridspec(nrows=6, ncols=4, left=0.05, right=0.48, wspace=0.15)


    proj_params = (ScoremapFOVProjParams & f'proj_param_id = "{proj_param_id}"').fetch1()

    for no, binning_param_set in enumerate(binning_parameters):

        print(f'Getting {binning_param_set} for animal {animal_name}')
        
        # Get entries and project / smooth
        score_maps, no_sessions, moran_is, scoremap_keys = get_scoremaps_moran(animal = animal_name, 
                                                                               binning_param_set = binning_param_set, 
                                                                               proj_param_id = proj_param_id,
                                                                               metasession_ref = metasession_ref,
                                                                               param_hash_id_cell = param_hash_id_cell,
                                                                               param_hash_session = param_hash_session,                                                       
                                                                               )
        print(f'Projection maps | statistic: {proj_params["proj_statistic"]} | sigma: {proj_params["proj_sigma"]}')
        projection, _, _, _, _ = stack_project_fov_score_maps(score_maps, 
                                                              statistic = proj_params['proj_statistic'], 
                                                              sigma = proj_params['proj_sigma'],
                                                              )


        # Prepare anatomical mask if needed 
        # and create overlay with retrieved projection
        if anat_mask_label is not None:
            print(f'Retrieving anatomical mask for region "{anat_mask_label}"')
            bins_x, bins_y = (ScoremapFOV & scoremap_keys).fetch1('bins_x','bins_y')
            scale_factor = np.mean([(bins_y[-1]-bins_y[0])/len(bins_y), (bins_x[-1]-bins_x[0])/len(bins_x)])
            bins_x, bins_y = bins_x.astype(int), bins_y.astype(int)
            
            # Retrieve anatomical masks 
            # The attributes 'mean_img', projection_short and 'center_plane' are just here to "sparsen" the output
            collected_masks_region = get_aligned_anatomical_masks(animal_name, 
                                                                  region=anat_mask_label, 
                                                                  projection_short='mean_img', 
                                                                  center_plane=0, 
                                                                  )
            collected_masks_region = collected_masks_region[bins_y[0]:bins_y[-1], bins_x[0]:bins_x[-1]]
            rescaled_masks         = rescale(collected_masks_region, 1/scale_factor)
            rescaled_masks         = rescaled_masks[:projection.shape[0], :projection.shape[1]] # Crop
            if anat_mask_sigma: 
                rescaled_masks     = gaussian_filter(rescaled_masks,  sigma=anat_mask_sigma)
            # Show anatomical mask
            projection = rescaled_masks #projection + (rescaled_masks / (rescaled_masks.max() / np.nanmax(projection)))





        # Prepare axes
        ax_scoremap = figure.add_subplot(gs[:-1, no])
        ax_moran    = figure.add_subplot(gs[-1:, no])
        

        # Draw Scoremap results
        cmap = getattr(binning_params, binning_param_set)['colormap']
        if anat_mask_label is None:
            vmin, vmax = getattr(binning_params, binning_param_set)['vmin_vmax']
        else:
            # Ignore if anatomical mask is to be plotted
            vmin, vmax = None, None

        g = ax_scoremap.imshow(projection, 
                               cmap=cmap, 
                               vmin=vmin, 
                               vmax=[np.nanpercentile(projection,percentile) if percentile is not None else vmax][0]
                               )
        
        
        # ... take care of colorbar
        divider = make_axes_locatable(ax_scoremap)
        cax = divider.append_axes("right", size="2%", pad=0.08)
        cbar = plt.colorbar(g, cax=cax)
        cbar.outline.set_visible(False)
        #cbar.ax.set_ylabel(f'{binning_param_set}', rotation=270, labelpad=14)
        #cbar.ax.get_yaxis().set_ticks([])
        ax_scoremap.get_xaxis().set_ticks([]); ax_scoremap.get_yaxis().set_ticks([])    
        ax_scoremap.set_title(f'{binning_param_set} | {no_sessions} sessions')
        # Draw Moran's I results
        try:
            ax_moran.hist(moran_is['moran_i_shuffles'], bins=50, lw=0, color='#333', alpha=1)
            ax_moran.axvline(x=moran_is['moran_i_95'], lw=1.5, color='k', ls=':', label='Moran\'s I 95th', alpha=.75)
            ax_moran.axvline(x=moran_is['moran_i_99'], lw=1.5, color='k', ls=':', label='Moran\'s I 99th')
            ax_moran.axvline(x=moran_is['moran_i'],    lw=2.5, color='r', ls='-', label='Moran\'s I Data')
        except ValueError:
            print('Morans I results can not be plotted')
        ax_moran.get_yaxis().set_ticks([])
    
    
    if bin_size_microns is not None:
        # Draw scalebar (defaults to 100. microns)
        scalebar_length_px = scalebar_length / bin_size_microns
        color_scalebar = [.2,.2,.2]
        xlim_scoremap, ylim_scoremap = ax_scoremap.get_xlim(), ax_scoremap.get_ylim()
        ax_scoremap.plot([np.max(xlim_scoremap)-scalebar_length_px-5,
                          np.max(xlim_scoremap)-scalebar_length_px-5+scalebar_length_px], 
                         [np.max(ylim_scoremap)-5,
                          np.max(ylim_scoremap)-5], 
                          lw=4, 
                          color=color_scalebar, 
                          alpha=.8, 
                          solid_capstyle='butt')

    sns.despine(left=True,bottom=False)
    plt.suptitle(f'{animal_name}', x=0.06, y=.95)

    if save_path is not None: 
        figure.savefig(save_path / f'scoremap collage {animal_name} {binning_parameters}.pdf', 
                       dpi=600, 
                       bbox_inches='tight')




def _create_additive_overlay(image1,image2):
    overlay = np.zeros((image1.shape[0], image1.shape[1],3))
    overlay[:,:,0] = (image2-image2.min())/(image2.max()-image2.min())
    overlay[:,:,1] = np.minimum(1, np.maximum(0,(image1-image2.min())/(image2.max()-image2.min())))
    overlay[:,:,2] = (image2-image2.min())/(image2.max()-image2.min())
    return overlay




def plot_summary_xcorr_scoremaps(animal_name,
                                 scoremap_entry_animal,
                                 invert_order=False,
                                 save_path=None,
                                 **kwargs,
                                ):
    '''
    Create single [!] summary plot
    for overlay of FOV score maps, and
    show results of for ScoremapCorr() (FOV scoremap cross correlations)

    
    Parameters
    ----------
    animal_name            :  str
                              animal_name
    scoremap_entry_animal  :  dict
                              Datajoint retrieval key for ScoremapCorr entry
    invert_order           :  boolean
                              Invert the order of scoremap A and B? 
    save_path              :  str or pathlib.Path
                              Save figure under this path
    
    **kwargs
        bin_size_microns : float
                           Bin size in microns. If given, scale bar will be displayed.
                           Defaults to None.
        scalebar_length  : float
                           Scale bar length in microns. CAVE will only be drawn
                           if bin_size_microns is not None.
                           Defaults to 100.


    '''
    bin_size_microns    = kwargs.get('bin_size_microns', None)
    scalebar_length     = kwargs.get('scalebar_length', 100.)

    if save_path is not None: 
        if isinstance(save_path, str):
            save_path = Path(save_path) 

    sns.set(style='white', font_scale=1.6)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    figure = plt.figure(constrained_layout=False, figsize=(10,11))
    gs = figure.add_gridspec(nrows=8, ncols=1, left=0.0, right=0.5, hspace=0.4)
   
    ### RETRIEVE PROJECTION PARAMETERS #####################################################################################
    statistic, sigma = (ScoremapFOVProjParams & scoremap_entry_animal).fetch1('proj_statistic', 'proj_sigma')


    ### RETRIEVE SCOREMAP ENTRIES ##########################################################################################
    score_maps_entry_A, score_maps_entry_B = get_scoremaps_a_b(scoremap_entry_animal)
    aligned_metas_A, aligned_metas_B       = score_maps_entry_A['aligned_metas'], score_maps_entry_B['aligned_metas']
    aligned_metas_A, aligned_metas_B       = prune_keys(aligned_metas_A, aligned_metas_B)

    ax        = figure.add_subplot(gs[:-1, 0])
    ax_xcorr  = figure.add_subplot(gs[-1:, 0])

    #### PLOT HISTOGRAM ####################################################################################################
    ax_xcorr.hist(scoremap_entry_animal['xcorr_shuffles'], bins=50, lw=0, color='#333', alpha=1);

    #### PLOT AXVLINE   ####################################################################################################
    ax_xcorr.axvline(x=scoremap_entry_animal['xcorr_95'], 
                     lw=1.5, color='k', 
                     ls=':', 
                     label='XCorr 95th', 
                     alpha=.75)
    ax_xcorr.axvline(x=scoremap_entry_animal['xcorr_99'], 
                     lw=1.5, color='k', 
                     ls=':', 
                     label='XCorr 99th')
    ax_xcorr.axvline(x=scoremap_entry_animal['xcorr_5'], 
                     lw=1.5, color='k', 
                     ls=':', 
                     label='XCorr 5th', 
                     alpha=.75)
    ax_xcorr.axvline(x=scoremap_entry_animal['xcorr_1'], 
                    lw=1.5, color='k', 
                    ls=':', 
                    label='XCorr 1st')
    ax_xcorr.axvline(x=scoremap_entry_animal['xcorr'], 
                     lw=1.5, color='red', 
                     ls=':', 
                     label='XCorr Data')
    ax_xcorr.get_yaxis().set_ticks([])


    #### CREATE PROJECTION ####################################################################################################

    projection_A, _,_,_,_ = stack_project_fov_score_maps(aligned_metas_A, statistic=statistic, sigma=sigma)
    projection_B, _,_,_,_ = stack_project_fov_score_maps(aligned_metas_B, statistic=statistic, sigma=sigma)
    
    if invert_order:
        projection_A_ = projection_A.copy()
        projection_A  = projection_B 
        projection_B  = projection_A_ 

    projection_A /= projection_A.max()
    projection_B /= projection_B.max()
    
    # OVERLAY 
    overlay = _create_additive_overlay(projection_A, projection_B)
    overlay = overlay[:,:,[1,0,2]].copy()

    ax.imshow(overlay)
    ax.get_xaxis().set_ticks([]);ax.get_yaxis().set_ticks([])

    if bin_size_microns is not None:
        # Draw scalebar (defaults to 100. microns)
        scalebar_length_px = scalebar_length / bin_size_microns
        color_scalebar = [1,1,1]
        xlim_scoremap, ylim_scoremap = ax.get_xlim(), ax.get_ylim()
        ax.plot([np.max(xlim_scoremap)-scalebar_length_px-5,
                 np.max(xlim_scoremap)-scalebar_length_px-5+scalebar_length_px], 
                [np.max(ylim_scoremap)-5,
                 np.max(ylim_scoremap)-5], 
                 lw=4, 
                 color=color_scalebar, 
                 alpha=.8, 
                 solid_capstyle='butt')


    # ... for title 
    binning_param_set_A = score_maps_entry_A["binning_param_set_A"]
    binning_param_set_B = score_maps_entry_B["binning_param_set_B"]
    n_A = len(aligned_metas_A)
    n_B = len(aligned_metas_B)
    ax.set_title(f'{binning_param_set_A} (n={n_A}) X {binning_param_set_B} (n={n_B})')
    
    sns.despine(left=True,bottom=False)
    if save_path is not None: 
        figure.savefig(save_path / f'scoremap collage {animal_name} {binning_param_set_A}x{binning_param_set_B}.pdf', 
                        dpi=600, bbox_inches='tight')





#######################################################################################################################################################
#######################################################################################################################################################
# Comparisons 

def plot_cumulative_horizontal(values_A, 
                               values_B, 
                               bins, 
                               labels_A=None, 
                               labels_B=None, 
                               color_A='cornflowerblue',
                               color_B='black',
                               xlabel='Cell size', 
                               Alabel='Grid animals', 
                               Blabel='OV animals', 
                               average=False,
                               median_mean='median',
                               save_path=None,
                               save_title=None):
    '''
    Create cumulative histogram plot and box plot from two arrays
    - values_A
    - values_B 
    and bins     
    
    labels: Len(values)
    '''
    assert median_mean in ['median','mean'], f'{median_mean} not valid. Choose either median or mean.'
    if median_mean=='median':
        meanmedian_A = np.nanmedian(values_A)
        meanmedian_B = np.nanmedian(values_B)
    else:
        meanmedian_A = np.nanmean(values_A)
        meanmedian_B = np.nanmean(values_B)

    bin_width = np.mean(np.diff(bins))
    center_bins = bins + bin_width/2
    center_bins = center_bins[:-1]
    
    if labels_A is not None: 
        assert len(labels_A) == len(values_A)
    if labels_B is not None:
        assert len(labels_B) == len(values_B)   
    if (labels_A is None) and (labels_B is None):
        labels_A = np.zeros_like(values_A)
        labels_B = np.zeros_like(values_B)
    
    # Initialize figure
    sns.set(style='white', font_scale=2.)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(3,4))
    gs = figure.add_gridspec(nrows=10, ncols=1, hspace=-.1)

    ax1 = figure.add_subplot(gs[:6, :])
    ax2 = figure.add_subplot(gs[7:, :])

    if average: 
        # A
        hist_values_A,_     = np.histogram(values_A, bins=bins)
        cdf_hist_values_A   = np.cumsum(hist_values_A)/np.cumsum(hist_values_A)[-1]
        ax1.plot(center_bins, cdf_hist_values_A, label=Alabel, color=color_A, alpha=.8, lw=3)
        # B
        hist_values_B,_     = np.histogram(values_B, bins=bins)
        cdf_hist_values_B   = np.cumsum(hist_values_B)/np.cumsum(hist_values_B)[-1]
        ax1.plot(center_bins, cdf_hist_values_B, label=Blabel, color=color_B, alpha=.7, lw=3)


    else: 
        for no,sub in enumerate(set(labels_A)):
            #label_A   = [f'{Alabel}' if not no else None][0]
            label_m_A = [f'{meanmedian_A:.1f}' if not no else None][0]
            
            hist_values_A,_     = np.histogram(values_A[labels_A==sub], bins=bins)
            cdf_hist_values_A   = np.cumsum(hist_values_A)/np.cumsum(hist_values_A)[-1]
            ax1.plot(center_bins, cdf_hist_values_A, label=label_m_A, color=color_A, alpha=.8, lw=3)


        for no,sub in enumerate(set(labels_B)):
            #label_B   = [f'{Blabel}' if not no else None][0]
            label_m_B = [f'{meanmedian_B:.1f}' if not no else None][0]
            
            hist_values_B,_     = np.histogram(values_B[labels_B==sub], bins=bins)
            cdf_hist_values_B   = np.cumsum(hist_values_B)/np.cumsum(hist_values_B)[-1]
            ax1.plot(center_bins, cdf_hist_values_B, label=label_m_B, color=color_B, alpha=.7, lw=3)
    
    # Take care of axis limits
    ax1.set_xlim(0,)
    xlims_ax1 = ax1.get_xlim()

    # Make *common* x tick labels
    xticks = np.arange(xlims_ax1[0],xlims_ax1[1],80)
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([])

    # Draw horizontal box plot underneath figure 

    boxprops_A    = dict(color=color_A,linewidth=1.5)
    medianprops_A = dict(color=color_A,linewidth=0)
    meanprops_A   = dict(color=color_A,linewidth=4.5, ls='-', solid_capstyle='butt')

    boxprops_B    = dict(color=color_B,linewidth=1.5)
    medianprops_B = dict(color=color_B,linewidth=0)
    meanprops_B   = dict(color=color_B,linewidth=4.5, ls='-', solid_capstyle='butt')

    # Plot box plot 
    ax2.boxplot(values_A, positions=[0],
                vert=False,
                widths=.8,
                showfliers=False,
                showmeans=True,
                meanline=True,
                whis=[1,99],
                boxprops=boxprops_A,
                medianprops=medianprops_A,
                meanprops=meanprops_A);
    ax2.boxplot(values_B, positions=[1],
                vert=False,
                widths=.8,
                showfliers=False,
                showmeans=True,
                meanline=True,
                whis=[1,99],
                boxprops=boxprops_B,
                medianprops=medianprops_B,
                meanprops=meanprops_B);

                
    ax2.set_yticks([0,1])
    ax2.set_yticklabels([Alabel, Blabel], rotation=0, ha='right')
    ax2.set_xlim(xlims_ax1)
    ax2.set_xticks(xticks)

    ax1.set_ylabel('Prob.')
    ax1.legend(loc='lower right')
    ax2.set_xlabel(xlabel)
   
    
    #plt.tight_layout()
    sns.despine(left=True)
    
    if save_path is not None:
        if isinstance(save_path, str):
            save_path = Path(save_path)
        print('Saving figure under {}'.format(str(save_path)))
        figure.savefig(save_path / f'{xlabel} {save_title}.pdf', bbox_inches='tight')
    
    #plt.show()
    


