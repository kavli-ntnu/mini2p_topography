### Figure 3C, D - Obenhaus et al. 
### Figure S3 A, B, and D - Obenhaus et al. 

import sys, os
import os.path
import numpy as np 
import datajoint as dj
import cmasher as cmr
from tabulate import tabulate
import pandas as pd 
from scipy.stats import pearsonr, spearmanr
from general import get_star

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from dj_schemas.utils import *
from dj_plotter.helpers.plotting_helpers import make_linear_colormap

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'

def _plot_summary(celltype_df, pure_only):
    '''
    FIGURE
    Create summary plot of cell type fraction (supplement).
    Plot the fraction of each major cell class against each other 
    in a grid of scatter plots

    '''
    sns.set(style='white',font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(10,9))

    ax = figure.add_subplot(3,3,1)
    # Grid OV 
    ax.scatter(celltype_df.ratio_grid, celltype_df.ratio_ovc, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('Grid'), ax.set_ylabel('OV')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    ax = figure.add_subplot(3,3,2)
    ax.scatter(celltype_df.ratio_grid, celltype_df.ratio_hd, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('Grid'), ax.set_ylabel('HD')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    ax = figure.add_subplot(3,3,3)
    ax.scatter(celltype_df.ratio_grid, celltype_df.ratio_border, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('Grid'), ax.set_ylabel('Border')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    # next row 
    ax = figure.add_subplot(3,3,4)
    ax.scatter(celltype_df.ratio_border, celltype_df.ratio_hd, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('Border'), ax.set_ylabel('HD')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    ax = figure.add_subplot(3,3,5)
    ax.scatter(celltype_df.ratio_ovc, celltype_df.ratio_hd, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('OV'), ax.set_ylabel('HD')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    ax = figure.add_subplot(3,3,6)
    ax.scatter(celltype_df.ratio_ovc, celltype_df.ratio_border, marker='o', s=90, color='k', lw=0, alpha=.6)
    ax.set_xlabel('OV'), ax.set_ylabel('Border')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    # # next row 
    # ax = figure.add_subplot(3,3,7)
    # ax.scatter(celltype_df.ratio_info, celltype_df.ratio_hd, marker='o', s=90, color='k', lw=0, alpha=.6)
    # ax.set_xlabel('Info'), ax.set_ylabel('HD')
    # ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    # ax = figure.add_subplot(3,3,8)
    # ax.scatter(celltype_df.ratio_info, celltype_df.ratio_grid, marker='o', s=90, color='k', lw=0, alpha=.6)
    # ax.set_xlabel('Info'), ax.set_ylabel('Grid')
    # ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    # ax = figure.add_subplot(3,3,9)
    # ax.scatter(celltype_df.ratio_info, celltype_df.ratio_ovc, marker='o', s=90, color='k', lw=0, alpha=.6)
    # ax.set_xlabel('Info'), ax.set_ylabel('OV')
    # ax.set_xlim(-.025,); ax.set_ylim(-.025,)




    plt.tight_layout()
    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'
    plt.savefig(
               f'{figure_dir}/fraction all cell types (pure={pure}).pdf', 
               dpi=300, 
               facecolor='w', 
               edgecolor='w',
               orientation='portrait',
               format=None,
               transparent=False, 
               bbox_inches='tight', 
               pad_inches=0,
              )
    #plt.show()


def _plot_grid_ov(celltype_df, pure_only, method='pearson', cmap='cmr.guppy'):
    '''
    FIGURE
    Create plot of cell type fractions for grid vs. ov cells.
    Plot the fraction of each major cell class against each other and 
    color code by animal name

    '''
    
    colors = make_linear_colormap(celltype_df.animal_name, categorical=True, cmap=cmap)
    

    # Only Grid and OV 
    ncelltypes_sub_nonan = celltype_df.copy().dropna()
    no_animals = len(set(ncelltypes_sub_nonan.animal_name))


    ratio_grid = ncelltypes_sub_nonan.ratio_grid
    ratio_ovc  = ncelltypes_sub_nonan.ratio_ovc 
    if method=='pearson':
        r, p = pearsonr(ratio_grid,ratio_ovc)
    elif method=='spearman':
        r, p = spearmanr(ratio_grid,ratio_ovc)
    else:
        raise NotImplementedError(f'Correlation method "{method}" not found')
        
    sig = get_star(p)

    # PRINT STATS
    print(f'Plotting grid vs. ov for {no_animals} animals across {len(ncelltypes_sub_nonan)} sessions')
    if method=='pearson':
        print(f'Pearsons r={r:.4} ngrid={len(ratio_grid)} novc={len(ratio_ovc)} p={p:.6f}{sig}')
    elif method=='spearman':
        print(f'Spearman r={r:.4} ngrid={len(ratio_grid)} novc={len(ratio_ovc)} p={p:.6f}{sig}')

    sns.set(style='white',font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(4,4))

    ax = figure.add_subplot(111)

    # Grid vs. OV    
    ax.scatter(celltype_df.ratio_grid, 
               celltype_df.ratio_ovc, 
               marker='o', 
               s=220, 
               c=colors, 
               lw=0, 
               alpha=.8,
               )

    xlim, ylim = ax.get_xlim(), ax.get_ylim()      
    ax.text(.7*xlim[1], .7*ylim[1], f'r: {r:.2f} {sig}') # Pearson's r
    ax.set_xlabel('Fraction Grid'), ax.set_ylabel('Fraction OV')
    ax.set_xlim(-.025,); ax.set_ylim(-.025,)

    divider = make_axes_locatable(ax)
    norm = mpl.colors.BoundaryNorm(np.arange(-.5, no_animals+.5, 1), eval(cmap).N) # Custom scalar mappable
    cax = divider.append_axes("right", size="5%", pad=0.08)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=eval(cmap)), 
                        cax=cax, 
                        ticks=np.arange(0,no_animals,1)
                        )
    
    cbar.outline.set_visible(False)

    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'

    plt.savefig(
            f'{figure_dir}/fraction grid vs. ov (pure={pure}) method={method}.pdf', 
            dpi=300, 
            facecolor='w', 
            edgecolor='w',
            orientation='portrait',
            format=None,
            transparent=False, 
            bbox_inches='tight', 
            pad_inches=0,
            )
    # plt.show()

def __return_p_pearsons(array_A,array_B):
    _, p = pearsonr(array_A, array_B)
    return p

def __return_p_spearmans(array_A,array_B):
    _, p = spearmanr(array_A, array_B)
    return p


def _plot_corr_heatmap(celltype_df, pure_only, method='pearson'):

    sns.set(style='white', font_scale=2.0)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    # Compute the correlation matrix
    corr = celltype_df.corr(method=method)

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))
    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(240, 28, s=100, l=70, as_cmap=True)

    figure = plt.figure(figsize=(5,5))
    ax = figure.add_subplot(111)
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-.7, vmax=.7, center=0,
                square=True, linewidths=.8, cbar_kws={"shrink": .5}, ax=ax)
    ax.set_yticks([0.5,1.5,2.5,3.5]);
    ax.set_yticklabels(['','OV','HD','Border'], rotation=0, ha='right', va='center');
    ax.set_xticks([0.5,1.5,2.5,3.5]);
    ax.set_xticklabels(['Grid','OV','HD','Border'], rotation=40, ha='right', va='top');

    ax.set_xlim(0,3)
    ax.set_ylim(4,0)

    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'

    plt.savefig(
            f'{figure_dir}/fraction correlations (pure={pure}) method={method}.pdf', 
            dpi=300, 
            facecolor='w', 
            edgecolor='w',
            orientation='portrait',
            format=None,
            transparent=False, 
            bbox_inches='tight', 
            pad_inches=0,
            )


    # Return view of corrs and p values 
    print('\nCorrelations')
    print(tabulate(corr, headers='keys', tablefmt='psql'))
    if method=='pearson':
        print('Pearson')
        p_values_corrs = celltype_df.corr(method=__return_p_pearsons)
    elif method=='spearman':
        print('Spearman')
        p_values_corrs = celltype_df.corr(method=__return_p_spearmans)
    print('\np values\n')
    print(tabulate(p_values_corrs, headers='keys', tablefmt='psql'))

    # plt.show()



def celltype_ratios(animals, 
                    pure_only=True,
                    param_hash_id_cell='standard',
                    param_hash_session='cf83e1357eefb8bd',
                    ov_cutoff_id='A',
                    region='MEC',
                    cutoff_n_cells=20
                    ):

    # Get cell parameter dictionary 
    cell_parameter_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')
    # Get rid of ov_cutoff_id -> set this custom 
    ov_cutoff_id_previous = cell_parameter_dict.pop('ov_cutoff_id')
    cell_parameter_dict['ov_cutoff_id'] = ov_cutoff_id
    print(f'Previous cutoff OV cells (ov_cutoff_id): {ov_cutoff_id_previous} | Set to: {ov_cutoff_id}\n')

    # Brain region filter 
    assert region in ['MEC','PAS'], f'Region "{region}" not understood. Choose "MEC" or "PAS"'
    if region == 'MEC':
        print('Choosing brain region filter MEC')
        brainregion = RoisCorrBrainLoc.MEC 
    elif region == 'PAS':
        print('Choosing brain region filter PAS')
        brainregion = RoisCorrBrainLoc.PAS

    # Filter cells by score / cutoff 
    bordercells = BorderScore * CutoffsBorderScore.proj('borderscore_95') & 'borderscore > borderscore_95'
    gridcells   = GridScore * CutoffsGridscore.proj('gridscore_95') & 'gridscore > gridscore_95'
    hdcells     = AngularRate.Stats * CutoffsMVL.proj('mvl_95') & 'mvl > mvl_95'
    #infocells   = Ratemap.Stats * CutoffsInfoContent.proj('info_content_95') & 'information_content > info_content_95'
    ovcells     = OVC.proj(is_ovc='is_ovc', session_name='base_session') & 'is_ovc > 0.5'
    
    # Build the cell type dictionary from scratch 
    # First just include all cells (do not filter for "pure")
    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )
    print(f'Found {len(all_sessions)} sessions\n')

    cell_no_dicts = []

    # ... for checking whether ov sessions exists or not 
    ovc_sessions = OVC.proj(is_ovc='is_ovc', session_name='base_session')

    # Loop over sessions and collect cell numbers / fractions
    animal_names  = []
    session_names = []
    for no, sess in enumerate(all_sessions.proj('animal_name')):

        filtered_cells = FilteredCells.proj() & sess & brainregion & f'param_hash_id_cell = "{param_hash_id_cell}"'
        if len(filtered_cells) < cutoff_n_cells:
            # Skip this session
            continue 
        else: 
            animal_names.append(sess['animal_name'])
            session_names.append(sess['session_name'])

        print(f'Processing {no+1}/{len(all_sessions)} (session_name {sess["session_name"]} - {len(filtered_cells)} cells)')
        
        border = bordercells & filtered_cells & cell_parameter_dict
        grid   = gridcells   & filtered_cells & cell_parameter_dict
        hd     = hdcells     & filtered_cells & cell_parameter_dict
        #info   = infocells   & filtered_cells & cell_parameter_dict
        ov     = ovcells     & filtered_cells & cell_parameter_dict
        
        if pure_only: 
            border_ = border.proj() - grid.proj() - hd.proj() - ov.proj()# - info.proj() 
            grid_   = grid.proj() - border.proj() - hd.proj() - ov.proj()# - info.proj() 
            hd_     = hd.proj() - border.proj() - grid.proj() - ov.proj()# - info.proj() 
            #info_   = info.proj() - hd.proj() - border.proj() - grid.proj() - ov.proj()
            ov_     = ov.proj() - border.proj() - grid.proj() - hd.proj()# - info.proj() 

            border  = border_ 
            grid    = grid_
            hd      = hd_
            #info    = info_
            ov      = ov_ 

        sess_dict = {}
        sess_dict['animal_name']  = sess['animal_name']
        sess_dict['session_name'] = sess['session_name']
        sess_dict['border']       = len(border)
        sess_dict['grid']         = len(grid)
        sess_dict['hd']           = len(hd)
        #sess_dict['info']         = len(info)
        sess_dict['total']        = len(filtered_cells)
        
        # Make sure you don't add zeros for non-existing object sessions
        if ovc_sessions & sess:
            sess_dict['ov']  = len(ov)
        else:
            sess_dict['ov']  = np.nan
            
        cell_no_dicts.append(sess_dict)
    
    # Print out quick statistics 
    animals_considered  = set(animal_names)
    sessions_considered = set(session_names)
    print(f'\n{len(animals_considered)} animals across {len(sessions_considered)} sessions\n')
    print(f'Animals: {list(animals_considered)}')
        
    # Create pandas dataframe from dict 
    cell_no_df = pd.DataFrame(cell_no_dicts)
    # Get ratios
    cell_no_df['ratio_grid']  = cell_no_df['grid']  / cell_no_df.total
    cell_no_df['ratio_ovc']   = cell_no_df['ov']    / cell_no_df.total
    cell_no_df['ratio_hd']    = cell_no_df['hd']    / cell_no_df.total
    #cell_no_df['ratio_info']  = cell_no_df['info']  / cell_no_df.total
    cell_no_df['ratio_border']= cell_no_df['border']/ cell_no_df.total


    ncelltypes_sub = cell_no_df[['animal_name','ratio_grid','ratio_ovc','ratio_hd','ratio_border']]

    print('animal_name','ratio_grid','ratio_ovc','ratio_hd','ratio_border')
    print(ncelltypes_sub.count())
    return ncelltypes_sub


if __name__ == "__main__":

    all_animals = [
               '90222','90218','90647',
               '82913','88592','89622',
               '87244','89841','60480',
               '87245','87187','88106',
               '94557','97045','97046',
               ] 


    print(f'Including {len(all_animals)} animals')

    # Cutoff number of cells 
    cutoff_n_cells = 15 
    ov_cutoff_id = 'A'

    # Retrieve cells, fractions and plot summary and grid vs. ov plot
    pure_only = False # Compare only "pure" cells to each other (no mixed selectivity)
    corr_method = 'spearman' # spearman

    if pure_only:
        print('Comparing only "pure" cells (no mixed selectivity)')
    else:
        print('Comparing cells above threshold (can have mixed selectivity)')

    ncelltypes_sub = celltype_ratios(all_animals, 
                                     pure_only=pure_only,
                                     param_hash_id_cell='standard',
                                     param_hash_session='cf83e1357eefb8bd',
                                     ov_cutoff_id = ov_cutoff_id,
                                     region='MEC',
                                     cutoff_n_cells=cutoff_n_cells)


    #### CREATE FIGURES ##############################################################################################
    print('Creating summary plot of cell class fractions')
    _plot_summary(ncelltypes_sub, pure_only)
    print('Creating plot of cell class fractions for grid vs. ov')
    _plot_grid_ov(ncelltypes_sub, pure_only, method=corr_method)
    print('Creating heat map of correlations')
    _plot_corr_heatmap(ncelltypes_sub, pure_only, method=corr_method)


    print(figure_dir)
    print('Success.')
  