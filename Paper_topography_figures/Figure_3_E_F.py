### Figure 3 E and F - Obenhaus et al. 

import sys, os
import os.path
import numpy as np 
import datajoint as dj
import pickle
import cmasher as cmr
from general import get_star, print_mannwhitneyu, print_wilcoxon

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import * 

##### GLOBAL PLOTTING OPTIONS ############################################################

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY/'


def _plot_gric_ov_frac(session_frac_gridcells, 
                       session_frac_ov_cells,
                       pure_only,
                       animal_names
                       ):
                       
    '''
    Plot cell type fractions (OV an Grid) - line plot 
    and show pie chart at the bottom, displaying the average.   
    
    
    '''
    sns.set(style='white', font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(constrained_layout=False, figsize=(2.1, 6))
    gs = figure.add_gridspec(nrows=7, ncols=1)

    #### LINE PLOT
    ax1 = figure.add_subplot(gs[:-3, 0])  
    
    for frac_ov, frac_grid in zip(session_frac_ov_cells, session_frac_gridcells):
        ax1.plot([1,2],[frac_ov, frac_grid], '-', color='k', lw=2, alpha=.8)
        ax1.scatter([1,2],[frac_ov, frac_grid], color='k', lw=0, s=100, alpha=.8)
        ax1.set_xticks([1,2])
        ax1.set_xticklabels(['OV','Grid'])
    
    ax1.set_ylabel('Fraction all cells')
    ax1.set_xlim(.8,2.2)
    ax1.set_ylim(0,)

    #### PIE CHART 
    ax2 = figure.add_subplot(gs[-2:, 0])
    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    sizes = [np.mean(session_frac_ov_cells), np.mean(session_frac_gridcells), 1-(np.mean(session_frac_ov_cells) + np.mean(session_frac_gridcells))]
    colors = ['#2BA9D9','#F66E5B','#ccc'] 
    labels = '', '', ''
    explode = (0.02,0.02,0.02)

    wedges, _, _ = ax2.pie(sizes, 
                           colors=colors, 
                           explode=explode, 
                           labels=labels, 
                           autopct='%1.1f%%', 
                           pctdistance=1.8,
                           shadow=False, 
                           startangle=0)
    for w in wedges:
        w.set_linewidth(0)
        w.set_edgecolor('k')

    # Draw circle
    centre_circle = plt.Circle((0,0), 0.50, fc='white')
    ax2.add_artist(centre_circle)
        
    ax2.axis('equal') 
    sns.despine(left=True,bottom=False)
    
    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'
    
    plt.savefig(
        f'{figure_dir}/fraction ov vs. grid (pure={pure}) {animal_names}.pdf', 
        dpi=300, 
        facecolor='w', 
        edgecolor='w',
        orientation='portrait',
        format=None,
        transparent=False, 
        bbox_inches='tight', 
        pad_inches=0,
        )
    return


def _plot_discrimination(session_disc,
                         pure_only,
                         animal_names
                         ):
    '''
    Boxplots of discrimination ratio of number of Grid and OV cells.
    (len(Grid)-len(OV)) / (len(Grid) + len(OV))

    Also calculates statistics (Wilcoxon) for difference against 0.
    '''
                       
    sns.set(style='white', font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    # Box plot properties 
    boxprops_B    = dict(color='#333', linewidth=0)
    meanprops_B   = dict(linewidth=0) 
    medianprops_B = dict(color='k', linewidth=2.5, ls='-', solid_capstyle='butt')
    whiskerprops_B= dict(color='#333', linewidth=1.5)
    capprops_B    = dict(color='#333')
    props_B = [boxprops_B, medianprops_B, meanprops_B, whiskerprops_B, capprops_B]

    figure = plt.figure(constrained_layout=False, figsize=(.7, 6))
    gs = figure.add_gridspec(nrows=7, ncols=1)

    # BOX PLOT
    ax1 = figure.add_subplot(gs[:-3, 0])  

    props = props_B
    bp = ax1.boxplot([session_disc],
                      patch_artist=True,
                      positions=[1], 
                      widths=.8,
                      showfliers=False,
                      whis=[1,99],
                      showmeans=False,
                      meanline=True,
                      boxprops=props[0],
                      medianprops=props[1],
                      meanprops=props[2],
                      whiskerprops=props[3],
                      capprops=props[4],
                      zorder=1,
                    )
    # Change color of box
    for box, color in zip(bp['boxes'],['#ccc']):
        box.set(facecolor=color)

    # Scatter
    for disc in session_disc:
        offset = (np.random.rand() - .5) / 3
        ax1.scatter(1 + offset, disc, s=70, color='k', alpha=.6, lw=0, zorder=2)

    ax1.set_xticks([1])
    ax1.set_xticklabels(['Disc. ratio'], rotation=90, ha='center')

    ax1.axhline(y=0,color='#333',ls=':',lw=1.5, alpha=.9)

    ax1.set_xlim(.5,1.5)
    ax1.set_ylim(-1.1,1.1)
    ax1.set_ylabel('(Grid-OV)/(Grid+OV)')

    sns.despine(left=True)

    # STATISTICS
    # Wilcoxon signed-rank test data vs. 0. 
    print_wilcoxon(np.array(session_disc) - 0., label='Disc.Index (against 0.)')

    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'
    
    plt.savefig(
        f'{figure_dir}/Disc ratio (pure={pure}) {animal_names}.pdf', 
        dpi=300, 
        facecolor='w', 
        edgecolor='w',
        orientation='portrait',
        format=None,
        transparent=False, 
        bbox_inches='tight', 
        pad_inches=0,
        )
    return


#### SHUFFLED ##############################################################################################################

def _plot_shuff_frac(session_frac_gridcells, 
                     session_frac_ov_cells,
                     shuffled_frac_animals,
                     pure_only,
                     animal_names
                     ):
    '''
    Boxplots of Data and shuffled results for the 
    fractions of OV and Grid cells 
       
    '''
                       
    sns.set(style='white', font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    # Box plot properties 
    boxprops_B    = dict(color='#333', linewidth=0)
    meanprops_B   = dict(linewidth=0) 
    medianprops_B = dict(color='k', linewidth=2.5, ls='-', solid_capstyle='butt')
    whiskerprops_B= dict(color='#333', linewidth=1.5)
    capprops_B    = dict(color='#333')
    props_B = [boxprops_B, medianprops_B, meanprops_B, whiskerprops_B, capprops_B]

    figure = plt.figure(constrained_layout=False, figsize=(1.8, 6))
    gs = figure.add_gridspec(nrows=7, ncols=1)

    #### BOX PLOT
    ax1 = figure.add_subplot(gs[:-3, 0])  

    props = props_B
    bp = ax1.boxplot([session_frac_ov_cells, 
                      shuffled_frac_animals['frac_ov_cells'],
                      session_frac_gridcells,
                      shuffled_frac_animals['frac_grid_cells'],                     
                     ],
            patch_artist=True,
            positions=[1,2,3.5,4.5], 
            widths=.8,
            showfliers=False,
            whis=[1,99],
            showmeans=False,
            meanline=True,
            boxprops=props[0],
            medianprops=props[1],
            meanprops=props[2],
            whiskerprops=props[3],
            capprops=props[4],
            )
    # Change color of boxes
    for box, label, color in zip(bp['boxes'], 
                                 ['ov','shuff_ov','grid','shuff_grid'], 
                                 ['#2BA9D9','#2BA9D9','#F66E5B', '#F66E5B'],
                                 ):
        box.set(facecolor=color)
        if 'shuff' in label: # Change transparency
            r_, g_, b_, _ = box.get_facecolor()
            box.set(facecolor=(r_, g_, b_, .3))    

    ax1.set_xticks([1,2,3.5,4.5])
    ax1.set_xticklabels(['OV','Shuffled OV','Grid','Shuffled Grid'], rotation=90, ha='center')

    ax1.set_xlim(.5,5.)
    ax1.set_ylim(0,)
    ax1.set_ylabel('Fraction all cells')

    sns.despine(left=True)

    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'
    
    plt.savefig(
        f'{figure_dir}/Shuffled frac ov vs. grid (pure={pure}) {animal_names}.pdf', 
        dpi=300, 
        facecolor='w', 
        edgecolor='w',
        orientation='portrait',
        format=None,
        transparent=False, 
        bbox_inches='tight', 
        pad_inches=0,
        )
    return


def _plot_shuff_discrimination(session_disc,
                               shuffled_frac_animals,
                               pure_only,
                               animal_names
                               ):
    '''
    Boxplots of Data and shuffled results for the 
    discrimination ratio of Grid and OV cells.
       
    '''
                       
    sns.set(style='white', font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    # Box plot properties 
    boxprops_B    = dict(color='#333', linewidth=0)
    meanprops_B   = dict(linewidth=0) 
    medianprops_B = dict(color='k', linewidth=2.5, ls='-', solid_capstyle='butt')
    whiskerprops_B= dict(color='#333', linewidth=1.5)
    capprops_B    = dict(color='#333')
    props_B = [boxprops_B, medianprops_B, meanprops_B, whiskerprops_B, capprops_B]

    figure = plt.figure(constrained_layout=False, figsize=(1.1, 6))
    gs = figure.add_gridspec(nrows=7, ncols=1)

    # BOX PLOT
    ax1 = figure.add_subplot(gs[:-3, 0])  

    props = props_B
    bp = ax1.boxplot([session_disc, 
                      shuffled_frac_animals['disc_grid_ov'],
                     ],
            patch_artist=True,
            positions=[1,2], 
            widths=.85,
            showfliers=False,
            whis=[1,99],
            showmeans=False,
            meanline=True,
            boxprops=props[0],
            medianprops=props[1],
            meanprops=props[2],
            whiskerprops=props[3],
            capprops=props[4],
            )
    # Change color of boxes
    for box, label, color in zip(bp['boxes'], 
                                 ['ratio','shuff_ratio'], 
                                 ['#ccc','#ccc'],
                                 ):
        box.set(facecolor=color)
        if 'shuff' in label: # Change transparency
            r_, g_, b_, _ = box.get_facecolor()
            box.set(facecolor=(r_, g_, b_, .3))    

    ax1.set_xticks([1,2])
    ax1.set_xticklabels(['Data','Shuffled'], rotation=90, ha='center')

    ax1.axhline(y=0, color='#333', ls=':', lw=1.5, alpha=.9)

    ax1.set_xlim(.5,2.5)
    ax1.set_ylim(-1.1,1.1)
    ax1.set_ylabel('(Grid-OV)/(Grid+OV)')

    sns.despine(left=True)

    if pure_only:
        pure = 'YES'
    else:
        pure = 'NO'
    
    plt.savefig(
        f'{figure_dir}/Shuffled disc ratio (pure={pure}) {animal_names}.pdf', 
        dpi=300, 
        facecolor='w', 
        edgecolor='w',
        orientation='portrait',
        format=None,
        transparent=False, 
        bbox_inches='tight', 
        pad_inches=0,
        )
    return



def grid_ov_frac(animals, 
                 shuffled_frac_animals,
                 pure_only=True,
                 param_hash_id_cell='standard',
                 param_hash_session='cf83e1357eefb8bd',
                 ov_cutoff_id='A',
                 region='MEC',
                 cutoff_n_cells=20):
    '''
    Create a line plot and pie chart underneath for cell type fractions
    per animal or multiple animals

    '''

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

    # Retrieve sessions
    all_sessions = (Session.proj('animal_name') * FilteredSessions 
                    & [f'animal_name = "{animal}"' for animal in animals]
                    & f'param_hash_session = "{param_hash_session}"'
                    )

    # All OVC entries
    ov_sessions = all_sessions & OVC.proj(session_name='base_session')


    session_frac_ov_cells =  []
    session_frac_gridcells = []
    session_disc = []

    animal_names  = []
    session_names = []
    for no, sess in enumerate(ov_sessions.proj('animal_name')):

        filtered_cells = FilteredCells.proj() & sess & brainregion & f'param_hash_id_cell = "{param_hash_id_cell}"'
        if len(filtered_cells) < cutoff_n_cells:
            # Skip this session
            continue 
        else: 
            animal_names.append(sess['animal_name'])
            session_names.append(sess['session_name'])
        
        print(f'Processing {no+1}/{len(ov_sessions)} (session_name {sess["session_name"]} - {len(filtered_cells)} cells)')

        # Filter further for Grid and OV
        ov_cells_session  = OVC.proj(session_name='base_session', is_ovc='is_ovc') & filtered_cells & cell_parameter_dict & 'is_ovc = 1'
        gridcells_session = GridScore * CutoffsGridscore & filtered_cells & cell_parameter_dict & 'gridscore > gridscore_95'
        
        if pure_only:
            # Clean up those mixed cells ... 
            ov_cells_session_   = ov_cells_session   - gridcells_session.proj()
            gridcells_session_  = gridcells_session  - ov_cells_session.proj()
            ov_cells_session    = ov_cells_session_
            gridcells_session   = gridcells_session_


        frac_ov_cells_session   = len(ov_cells_session)  / len(filtered_cells)
        frac_grid_cells_session = len(gridcells_session) / len(filtered_cells)

        session_frac_ov_cells.append(frac_ov_cells_session)
        session_frac_gridcells.append(frac_grid_cells_session)

        # Discrimination index
        disc = (len(gridcells_session)-len(ov_cells_session)) / (len(gridcells_session)+len(ov_cells_session))
        session_disc.append(disc)

        print(f'{sess["session_name"]} | Frac OV: {frac_ov_cells_session:.2f} | Frac Grid: {frac_grid_cells_session:.2f} | Disc. index: {disc:.2f}')
    
    # Quick stats for display: 
    animals_considered  = set(animal_names)
    sessions_considered = set(session_names)
    print(f'\n{len(animals_considered)} animals across {len(sessions_considered)} sessions\n')
    print(f'Animals: {list(animals_considered)}')
    
    # Process shuffled results 
    if shuffled_frac_animals is not None:
        shuffled_frac_grid  = shuffled_frac_animals['frac_grid_cells']
        shuffled_frac_ov    = shuffled_frac_animals['frac_ov_cells']
        shuffled_disc       = shuffled_frac_animals['disc_grid_ov']
    else: 
        shuffled_frac_grid, shuffled_frac_ov, shuffled_disc = None, None, None
    
    # Print statistics 
    _get_stats(session_frac_gridcells, 
               session_frac_ov_cells, 
               session_disc,
               shuffled_frac_grid, 
               shuffled_frac_ov,
               shuffled_disc)

    return session_frac_gridcells, session_frac_ov_cells, session_disc


def _get_stats(session_frac_gridcells, 
               session_frac_ov_cells, 
               session_disc,
               shuffled_frac_grid, 
               shuffled_frac_ov,
               shuffled_disc
               ):

    # Report some more stats on fraction grid / OV cells and their ratios
    mean_frac_grid, sd_frac_grid = np.nanmean(session_frac_gridcells), np.nanstd(session_frac_gridcells)
    mean_frac_ov, sd_frac_ov     = np.nanmean(session_frac_ov_cells), np.nanstd(session_frac_ov_cells)
    mean_disc, sd_disc           = np.nanmean(session_disc), np.nanstd(session_disc)

    print_mannwhitneyu(session_frac_gridcells, 
                       session_frac_ov_cells,
                       'Grid',
                       'OV')

    print(f'\nGrid mean ± SD:         {mean_frac_grid:.4f} ± {sd_frac_grid:.4f}')
    print(f'OV mean ± SD:           {mean_frac_ov:.4f} ± {sd_frac_ov:.4f}')
    print(f'Disc. ratio mean ± SD:  {mean_disc:.4f} ± {sd_disc:.4f}\n')


    if (shuffled_frac_grid is not None) and (shuffled_frac_ov is not None) and (shuffled_disc is not None): 
        # Report some more stats on shuffled fraction grid / OV cells
        shuff_mean_frac_grid, shuff_sd_frac_grid = np.nanmean(shuffled_frac_grid), np.nanstd(shuffled_frac_grid)
        shuff_mean_frac_ov, shuff_sd_frac_ov     = np.nanmean(shuffled_frac_ov),  np.nanstd(shuffled_frac_ov)
        shuff_mean_disc, shuff_sd_disc           = np.nanmean(shuffled_disc),  np.nanstd(shuffled_disc)
        print(f'SHUFFLED Grid mean ± SD:           {shuff_mean_frac_grid:.4f} ± {shuff_sd_frac_grid:.4f}')
        print(f'SHUFFLED OV mean ± SD:             {shuff_mean_frac_ov:.4f} ± {shuff_sd_frac_ov:.4f}')
        print(f'SHUFFLED Disc. ratio mean ± SD:    {shuff_mean_disc:.4f} ± {shuff_sd_disc:.4f}\n')

        # Wilcoxon signed-rank test data vs. shuffled
        print('-- Wilcoxon signed-rank test data vs. mean of shuffled distribution:')
        print_wilcoxon(session_frac_gridcells - shuff_mean_frac_grid, 'Grid')
        print_wilcoxon(session_frac_ov_cells - shuff_mean_frac_ov, 'OV')
        # Compare Disc.Index
        print_wilcoxon(session_disc - shuff_mean_disc, 'Disc.Index')


if __name__ == "__main__":
 
    grid_mice = [ 
                '82913','88592', '87244', '60480',
                '97046','89841' 
                ]
    ov_mice   = [ 
                '87187','88106','87245','90222',
                '94557','89622'
                ]

    animals = ov_mice # grid_mice # ['89622']  #ov_mice #['85590'] #grid_mice # '82913' 

    print(f'Creating OV vs. Grid figure for {len(animals)} animal(s)')
    print(animals)
    print('\n')

    # Cutoff number of cells 
    cutoff_n_cells = 15 
    ov_cutoff_id = 'A'
    pure_only = True # Compare only "pure" cells to each other (no mixed selectivity)

    if pure_only:
        print('Comparing only "pure" cells (no mixed selectivity)')
    else:
        print('Comparing cells above threshold (can have mixed selectivity')


    # Shuffled results 
    # Load only if mice are either OV or GRID mice 
    # This is in `fig_2_E_F_Shuffling.py`
    dict_export = 'data_exports/shuffled_fig_3_F'
    # Load dictionaries with keys 
    # - frac_ov_cells
    # - frac_grid_cells
    # - disc_grid_ov # i.e. discrimination ratio (G-OV) / (G+OV)
    if animals == ov_mice: 
        with open(dict_export + f'/shuffled_ov_animals_frac_pure={pure_only}ov_cutoff_id={ov_cutoff_id}.pkl', 'rb') as f:
            shuffled_frac_animals =  pickle.load(f)
    elif animals == grid_mice: 
        with open(dict_export + f'/shuffled_grid_animals_frac_pure={pure_only}ov_cutoff_id={ov_cutoff_id}.pkl', 'rb') as f:
            shuffled_frac_animals =  pickle.load(f)
    else:
        shuffled_frac_animals = None
        print('No shuffled results were loaded.')

    session_frac_gridcells, session_frac_ov_cells, session_disc = grid_ov_frac(animals, 
                                                                               shuffled_frac_animals,
                                                                               pure_only=pure_only,
                                                                               param_hash_id_cell='standard',
                                                                               param_hash_session='cf83e1357eefb8bd',
                                                                               ov_cutoff_id = ov_cutoff_id,
                                                                               region='MEC',
                                                                               cutoff_n_cells=cutoff_n_cells)

    # Export the stuff - want to play with it
    with open(dict_export + f'/session_frac_gridcells_pure={pure_only}{animals}.pkl', 'wb') as f:
        pickle.dump(session_frac_gridcells, f, pickle.HIGHEST_PROTOCOL)
    with open(dict_export + f'/session_frac_ov_cells={pure_only}{animals}.pkl', 'wb') as f:
        pickle.dump(session_frac_ov_cells, f, pickle.HIGHEST_PROTOCOL)
    with open(dict_export + f'/session_disc={pure_only}{animals}.pkl', 'wb') as f:
        pickle.dump(session_disc, f, pickle.HIGHEST_PROTOCOL)


    #### CREATE FIGURES ##############################################################################################################################
    # LINE PLOT and PIE CHART
    _plot_gric_ov_frac(session_frac_gridcells, 
                       session_frac_ov_cells, 
                       pure_only, 
                       animals,
                       )
    # BOX PLOTS
    _plot_discrimination(session_disc,
                         pure_only,
                         animals
                         )


    if shuffled_frac_animals is not None: 
        # Plot box plots for data vs. shuffled 
        # Fractions: Comparison against shuffled distribution
        _plot_shuff_frac(session_frac_gridcells, 
                         session_frac_ov_cells, 
                         shuffled_frac_animals,  
                         pure_only, 
                         animals,
                         )
        # Disc. ratios: Comparison against shuffled distribution
        _plot_shuff_discrimination(session_disc,
                                   shuffled_frac_animals,
                                   pure_only,
                                   animals,
                                   )
    print(figure_dir)
    print('Success.')