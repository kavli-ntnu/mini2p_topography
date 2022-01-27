### Figure 3 Obenhaus et al. 
# 
# Panels:
# A and B 
#
# Cell maps and example ratemaps 

import sys, os
import os.path
import numpy as np 
import datajoint as dj
import cmasher as cmr

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from dj_plotter.helpers.plotting_helpers import *
from dj_plotter import dj_plotter

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### GLOBAL PLOTTING OPTIONS ############################################################
dot_size_rois = 2.5
scale_bar_rois = 50
draw_numbers = False

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'

def plot_cellmaps(session_name,
                  param_hash_id_cell='standard',
                  color_non_grid = [.7,.7,.7],
                  color_non_hd =   [.7,.7,.7],
                  color_non_border=[.7,.7,.7],
                  color_non_bvs=   [.7,.7,.7],
                  color_non_ovc=   [.7,.7,.7],
                  plot_ov = True,
                  **kwargs,
                 ):
    '''
    Plot cell maps as on display in Figure 3A.
    Cell (score) maps for Grid and OV 

    '''
    draw_pixels  = kwargs.get('draw_pixels',True)
    draw_centers = kwargs.get('draw_centers',False)
    dot_size     = kwargs.get('dot_size', dot_size_rois)

    filtered_cells = FilteredCells & f'session_name = "{session_name}"'  & f'param_hash_id_cell = "{param_hash_id_cell}"'
    cells = Session * ProjectionCorr * RoisCorr & filtered_cells 

    # Get cell parameter dictionary 
    cell_parameter_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')

    print(f'\nFound {len(cells)} cells for session {session_name}\n')


    gridscore, gridscore_95               = (GridScore * CutoffsGridscore 
                                             & cells.proj() 
                                             & cell_parameter_dict).fetch('gridscore','gridscore_95')

    mvl, mvl_95                           = (AngularRate.Stats * CutoffsMVL 
                                             & cells.proj() 
                                             & cell_parameter_dict).fetch('mvl','mvl_95')

    borderscore, borderscore_95           = (BorderScore * CutoffsBorderScore 
                                             & cells.proj() 
                                             & cell_parameter_dict).fetch('borderscore','borderscore_95')

    bvs, bvs_cutoff                       = (BVScore 
                                             & cells.proj() 
                                             & cell_parameter_dict).fetch('bvs'), .56

    if plot_ov:
        ovscore, ovscore_95, is_ovc       = (OVC * CutoffsOVScore 
                                             & cells.proj(base_session='session_name') 
                                             & cell_parameter_dict).fetch('ovscore','ovscore_95','is_ovc')
        ovscore = np.nan_to_num(ovscore)


    gridscore = np.nan_to_num(gridscore)
    mvl = np.nan_to_num(mvl)
    borderscore = np.nan_to_num(borderscore)
    bvs = np.nan_to_num(bvs)
    
    # Quick stats
    print(f'Grids  n={np.sum([gridscore>gridscore_95])}')
    print(f'HD     n={np.sum([mvl>mvl_95])}')
    print(f'Border n={np.sum([borderscore>borderscore_95])}')
    print(f'BVS    n={np.sum([bvs>bvs_cutoff])}')
    if plot_ov:
        print(f'OV     n={np.sum([is_ovc>.5])}')


    ##### DEFINE COLORBARS #################################################################################################

    gs_min, gs_max = 0, 2
    mvl_min, mvl_max = 0, .8
    bord_min, bord_max = 0, 1.5
    bvs_min, bvs_max = 0, 1
    if plot_ov:
        _, ovscore_max = 0, 1

    cmap_grid = 'Blues' #'Reds' #alternative: cmr.ember.colors    
    cmap_mvl  = cmr.cosmic.colors   
    cmap_bvs  = cmr.lavender.colors 
    cmap_ovc  = 'Blues' # alternative: cmr.cosmic.colors   

    colors_gs = make_linear_colormap(gridscore, 
                                     cmap=cmap_grid, 
                                     reference_numbers=np.array([gs_min, gs_max]))
    colors_gs[gridscore < gridscore_95] = color_non_grid

    colors_mvl = make_linear_colormap(mvl, 
                                     cmap=cmap_mvl, 
                                     reference_numbers=np.array([mvl_min, mvl_max]))
    colors_mvl[mvl < mvl_95] = color_non_hd

    colors_border = make_linear_colormap(borderscore, 
                                         cmap=cmap_bvs, 
                                         reference_numbers=np.array([bord_min, bord_max]))
    colors_border[borderscore < borderscore_95] = color_non_border

    colors_bvs = make_linear_colormap(bvs, 
                                      cmap=cmap_bvs, 
                                      reference_numbers=np.array([bvs_min, bvs_max]))
    colors_bvs[bvs < bvs_cutoff] = color_non_bvs
    if plot_ov: 
        colors_ov = make_linear_colormap(ovscore, 
                                        cmap=cmap_ovc, 
                                        reference_numbers=np.array([0, ovscore.max()]))
        colors_ov[is_ovc < .5] = color_non_ovc # < .5 because `is_ovc` is either 0 or 1 

    ##### INITIALISE PLOTTING ###############################################################################################

    plot = dj_plotter(cells, save_path=figure_dir, save_format='png')



    ##### GRIDS #############################################################################################################
    print('Plotting grids ...')
    plot.rois(
              draw_image=False, 
              draw_pixels=draw_pixels, 
              draw_numbers=draw_numbers,
              draw_centers=draw_centers, 
              draw_outlines=False, 
              colors=colors_gs, 
              display_title=False,
              dot_size=dot_size, 
              scalebar=scale_bar_rois, 
              path_suffix=f'{session_name} grid'
              )

    # Create a matching colorbar where values that are below threshold are greyed out
    cmap_ = sns.color_palette(cmap_grid, 256, 1)
    idx_grey = int((gridscore_95[0]/gs_max) * len(cmap_))
    colors_cmap = list(cmap_)
    colors_cmap_thresh = []

    for no, color in enumerate(colors_cmap):
        if no > idx_grey:
            colors_cmap_thresh.append(color)
        else:
            colors_cmap_thresh.append([.6, .6, .6])

    # Save colorbar for plot
    cb = make_colorbar(np.arange(1),
                  no_steps=len(colors_cmap_thresh), 
                  cmap=colors_cmap_thresh,
                  return_figure=True)
    cb.savefig(
               f'{figure_dir}/colorbar grid {session_name}.png', 
               dpi=300, 
               facecolor='w', 
               edgecolor='w',
               orientation='portrait',
               format=None,
               transparent=False, 
               bbox_inches='tight', 
               pad_inches=0,
              )



    if plot_ov:

        ##### OV ##############################################################################################################
        print('Plotting ovs ...')
        plot.rois(
                draw_image=False, 
                draw_pixels=draw_pixels, 
                draw_numbers=draw_numbers,
                draw_centers=draw_centers, 
                draw_outlines=False, 
                colors=colors_ov, 
                display_title=False,
                dot_size=dot_size, 
                scalebar=scale_bar_rois, 
                path_suffix=f'{session_name} ov'
                )


        # Create a matching colorbar where values that are below threshold are greyed out
        cmap_ = sns.color_palette(cmap_ovc, 256, 1)
        idx_grey = int((ovscore_95[0]/ovscore_max) * len(cmap_))
        colors_cmap = list(cmap_)
        colors_cmap_thresh = []

        for no, color in enumerate(colors_cmap):
            if no > idx_grey:
                colors_cmap_thresh.append(color)
            else:
                colors_cmap_thresh.append([.6, .6, .6])

        # Save colorbar for plot
        cb = make_colorbar(np.arange(1),
                    no_steps=len(colors_cmap_thresh), 
                    cmap=colors_cmap_thresh,
                    return_figure=True)
        cb.savefig(
                f'{figure_dir}/colorbar ov {session_name}.png', 
                dpi=300, 
                facecolor='w', 
                edgecolor='w',
                orientation='portrait',
                format=None,
                transparent=False, 
                bbox_inches='tight', 
                pad_inches=0,
                )

    print(f'\nFigures exported successfully.')
    return








if __name__ == "__main__":

    grid_session_name = "59825ec5641c94b4" # Grid animal - 88592
    ov_session_name   = "0ed19ecd643fdafa" # OV animal - 87245
    param_hash_id_cell = 'standard' # FilteredCells ! 

    # Cell maps
    plot_cellmaps(grid_session_name, param_hash_id_cell=param_hash_id_cell, color_non_ovc=[1.,1.,1.])
    plot_cellmaps(ov_session_name, param_hash_id_cell=param_hash_id_cell, color_non_ovc=[1.,1.,1.])



    # Plot example ratemaps
    cell_parameter_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')

    # Grid 
    print('\nExporting Grid examples')
    grid_cell_examples = [237,124,123,145] # 59825ec5641c94b4
    rm_entries         = (Session.proj() * Ratemap 
                       & f'session_name = "{grid_session_name}"' 
                       & cell_parameter_dict 
                       & [f'cell_id = "{id_}"' for id_ in grid_cell_examples])

    plot = dj_plotter(rm_entries, 
                      plots_per_view=1, 
                      save_path=figure_dir, 
                      save_format='png')
                      
    plot.ratemaps(hash_or_animal='hash', cue_card_pos='west', display_title=True, cmap='cmr.torch')

    # OV
    # This needs partial connection to the Moser storage to plot correctly. 
    # However, the information is available as is in the database (no additional storage conn. required)
    #print('\nExporting OV examples')
    # ov_cell_examples = [1,2,86,8] # 0ed19ecd643fdafa
    # rm_entries       = (Session.proj(base_session='session_name') * Ratemap.proj(base_session='session_name') 
    #                  * (OVC & 'is_ovc = 1').proj(ovscore='ovscore')
    #                  & f'base_session = "{ov_session_name}"' 
    #                  & cell_parameter_dict 
    #                  & [f'cell_id = "{id_}"' for id_ in ov_cell_examples])

    # plot = dj_plotter(rm_entries,
    #                   save_path=figure_dir,
    #                   save_format='pdf')

    # plot.ratemaps_ov(hash_or_animal='hash', display_title=True, cmap='cmr.torch')
    # print(figure_dir)
    # print('Success.')