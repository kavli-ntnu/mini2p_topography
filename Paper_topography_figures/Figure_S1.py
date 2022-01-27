### Figure S1 Obenhaus et al. 

# Panels: 
# B: FOV sizes
# C: SNR ratios 
# D: Number of cells


import sys, os
import os.path
from pathlib import Path
import numpy as np 
import datajoint as dj

# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# Prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

##### IMPORTS ###########################################################################
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import pandas as pd 
from scipy.stats import mannwhitneyu
from general import get_star

##### LOAD SCHEMA COMPONENTS #############################################################
from dj_schemas.dj_conn import *

##### EXPORT LOCATION ####################################################################
figure_dir = 'YOUR_EXPORT_DIRECTORY'
figure_dir = Path(figure_dir)


def _get_cell_count(animals, 
                    animals_trans,
                    animals_virus,
                    param_hash_id_cell = 'standard',
                    param_hash_session = 'cf83e1357eefb8bd',
                    min_no_cells = 10,
                    ):
    '''
    Get number of FilteredCells() per session and plot

    '''

    cells_  = (Session.proj('animal_name') * FilteredCells.proj() 
                & FilteredSessions 
                & f'param_hash_id_cell = "{param_hash_id_cell}"'
                & f'param_hash_session = "{param_hash_session}"'
                & [f'animal_name = "{animal}"' for animal in animals]
                )
    cells_df = pd.DataFrame(cells_.fetch(as_dict=True))

    label_trans_virus = []
    for ani in cells_df['animal_name'].values:
        if ani in animals_trans:
            label_trans_virus.append('Transg')
        elif ani in animals_virus: 
            label_trans_virus.append('Virus')
        else:
            label_trans_virus.append(None)
    cells_df['label'] = label_trans_virus

    cells_df = cells_df.groupby(by=['session_name','label']).count()
    cells_df.reset_index(inplace=True)
    cells_df = cells_df[['session_name','label','cell_id']]                
    # Filter by minimum number of cells 
    cells_df = cells_df.loc[cells_df.cell_id>min_no_cells]

    virus_avg = np.nanmean(cells_df[cells_df.label == 'Virus'].cell_id.values)
    virus_std = np.nanstd(cells_df[cells_df.label == 'Virus'].cell_id.values)
    trans_avg = np.nanmean(cells_df[cells_df.label == 'Transg'].cell_id.values)
    trans_std = np.nanstd(cells_df[cells_df.label == 'Transg'].cell_id.values)
    print(f'# of sessions: {len(cells_df)}')
    print(f'Virus average: {virus_avg:.2f} | SD: {virus_std:.2f}')
    print(f'Transg average: {trans_avg:.2f} | SD: {trans_std:.2f}')

    trans = cells_df[cells_df.label == 'Transg'].cell_id.values 
    virus = cells_df[cells_df.label == 'Virus'].cell_id.values
    mw_u, mw_p = mannwhitneyu(trans, virus)
    sig = get_star(mw_p)
    print(f'Mann-Whitney U={mw_u:.3f}, nTransg={len(trans)}, nVirus={len(virus)}, p={mw_p:.6f}{sig} two-sided')

    sns.set(style='white',font_scale=1.8)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(figsize=(4,2))
    pal = {"Transg": "#549287", "Virus": "#54D387"}
    sns.boxplot(data=cells_df,
                x='cell_id',
                y='label', boxprops=dict(alpha=1, linewidth=3), 
                width=0.5, 
                linewidth=2, 
                color='w', 
                palette=pal,
                orient='h')

    p1 = plt.axvline(x=virus_avg, color='#54D387', label=f'Virus {virus_avg:.1f}')
    p2 = plt.axvline(x=trans_avg, color='#549287', label=f'Transg {trans_avg:.1f}')
    sns.despine(left=True,bottom=False)

    plt.ylabel('')
    plt.xlabel('Cells/session')

    plt.legend(handles=[p1,p2], title='', bbox_to_anchor=(1.05, 1), loc='upper left')
    figure.savefig(figure_dir / f'figure cell numbers.pdf', bbox_inches='tight')
    return


def _get_snr_virus_trans(animals_trans,
                         animals_virus,
                         param_hash_id_cell ='standard',
                         param_hash_session ='cf83e1357eefb8bd',
                         ):
    '''
    Compare SNR (signal to noise ratio) between virus injected
    and transgenic animals

    '''
    animals = animals_trans + animals_virus
    print(f'{len(animals_trans)} transgenic animals and')
    print(f'{len(animals_virus)} virus injected animals')

    # Get cell entries filter ... FilteredCellsParams
    cell_params_dict = (FilteredCellsParams & f'param_hash_id_cell = "{param_hash_id_cell}"').fetch1('parameter_dict_cell')


    # DO NOT filter by FilteredCells() here, since those are after SNR filtering (you want the complete distribution)
    snr_df_f = (Session.proj('animal_name') * SNR.proj('snr_df_f') 
                    & (FilteredSessions.proj() & f'param_hash_session = "{param_hash_session}"')
                    & cell_params_dict
                    & [f'animal_name = "{animal}"' for animal in animals]
                    )
    cells_df = pd.DataFrame(snr_df_f.fetch(as_dict=True))


    label_trans_virus = []
    for ani in cells_df['animal_name'].values:
        if ani in animals_trans:
            label_trans_virus.append('Transg')
        elif ani in animals_virus: 
            label_trans_virus.append('Virus')
        else:
            label_trans_virus.append(None)
    cells_df['label'] = label_trans_virus
 
    snr_df_f_trans = cells_df[cells_df.label == 'Transg'].snr_df_f.values
    snr_df_f_virus = cells_df[cells_df.label == 'Virus'].snr_df_f.values

    # Prepare distribution plots
    bins = np.linspace(0,15,60)
    center_bins = bins + np.mean(np.diff(bins))/2
    center_bins = center_bins[:-1]

    hist_snr_virus, _  = np.histogram(snr_df_f_virus, bins=bins)
    hist_snr_trans, _  = np.histogram(snr_df_f_trans, bins=bins)
    cdf_hist_snr_virus = np.cumsum(hist_snr_virus)/np.cumsum(hist_snr_virus)[-1]
    cdf_hist_snr_trans = np.cumsum(hist_snr_trans)/np.cumsum(hist_snr_trans)[-1]
        
    # PLOT DISTRIBUTION
    sns.set(style='white',font_scale=1.8)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    figure = plt.figure(figsize=(7,6))
    gs = figure.add_gridspec(nrows=2, ncols=12, wspace=200, hspace=.5)


    ax = figure.add_subplot(gs[0,:6])
    ax.plot(center_bins, cdf_hist_snr_virus, color='#54D387', alpha=.8,lw=2)
    ax.plot(center_bins, cdf_hist_snr_trans, color='#549287', alpha=.8,lw=2)

    ax.set_ylabel('Like. of \noccurence')
    ax.set_xlabel('SNR')

    ax = figure.add_subplot(gs[0,7:])
    ax.step(center_bins, hist_snr_virus, color='#54D387', alpha=.7, lw=2)
    ax.step(center_bins, hist_snr_trans, color='#549287', alpha=.7, lw=2)
    ax.set_ylabel('Count')
    ax.set_xlabel('SNR')

    median_virus = np.nanmedian(snr_df_f_virus)
    median_trans = np.nanmedian(snr_df_f_trans)
    ax.axvline(x=3.5, color='#444', label='Cutoff', ls='--')
    ax.axvline(x=median_virus, color='#54D387', label='AAV', ls='--')
    ax.axvline(x=median_trans, color='#549287', label='Transg.', ls='--')
    
    print(f'# of cells: Virus: {len(snr_df_f_virus)} | Transg: {len(snr_df_f_trans)}')
    print(f'SNR virus  animals (median): {median_virus:.2f}')
    print(f'SNR transg animals (median): {median_trans:.2f}')
    
    sns.despine(left=True)


    # PLOT SESSION AVERAGE
    ax = figure.add_subplot(gs[1,:7])
    cells_df = cells_df.groupby(by=['session_name','label']).median()
    print(f'Session median over {len(cells_df)} sessions')
    cells_df.reset_index(inplace=True)
    cells_df = cells_df[['session_name','label','snr_df_f']]        


    pal = {"Transg": "#549287", "Virus": "#54D387"}
    sns.boxplot(data=cells_df,
                x='snr_df_f',
                y='label',boxprops=dict(alpha=1, linewidth=3), 
                width=0.5, 
                linewidth=2, 
                color='w', 
                palette=pal,
                orient='h',
                ax=ax)

    p1 = plt.axvline(x=median_virus, color='#54D387', label=f'Virus {median_virus:.1f}', ls='--')
    p2 = plt.axvline(x=median_trans, color='#549287', label=f'Transg {median_trans:.1f}', ls='--')
    plt.legend(handles=[p1,p2], title='', bbox_to_anchor=(1.05, 1), loc='upper left')

    sns.despine(left=True, bottom=False)

    plt.ylabel('')
    plt.xlabel('SNR / session')

    figure.savefig(figure_dir / 'SNR transgenic vs. virus session.pdf', bbox_inches='tight')

    # STATS
    trans = cells_df[cells_df.label == 'Transg'].snr_df_f.values
    virus = cells_df[cells_df.label == 'Virus'].snr_df_f.values
    mw_u, mw_p = mannwhitneyu(trans, virus, alternative='two-sided')
    sig = get_star(mw_p)
    print(f'Mann-Whitney U={mw_u:.3f}, nTransg={len(trans)}, nVirus={len(virus)}, p={mw_p:.6f}{sig} two-sided')

    return



def _get_fov_size(animals,
                  param_hash_session ='cf83e1357eefb8bd'
                  ):
    '''
    Plot FOV sizes, 
    which are the "width_microns_eff" and "height_microns_eff" entries 
    in ProjectionCorr

    '''

    sessions = (Session.proj('animal_name')
                    & (FilteredSessions.proj() & f'param_hash_session = "{param_hash_session}"')
                    & [f'animal_name = "{animal}"' for animal in animals]
                    )

    min_width = min_height = 250
    width_microns_eff, height_microns_eff = (ProjectionCorr 
                                                & sessions
                                                & f'width_microns_eff  > {min_width}'
                                                & f'height_microns_eff > {min_height}'
                                                & 'center_plane=0'
                                                ).fetch('width_microns_eff','height_microns_eff')


    print(f'FOV sizes over {len(width_microns_eff)} sessions')
    mean_width  = np.nanmean(width_microns_eff)
    mean_height = np.nanmean(height_microns_eff)
    std_width  = np.nanstd(width_microns_eff)
    std_height = np.nanstd(height_microns_eff)
    print(f'Mean width  (microns) ± SD: {mean_width:.2f} ± {std_width:.2f}')
    print(f'Mean height (microns) ± SD: {mean_height:.2f} ± {std_height:.2f}')


    # PLOT
    sns.set(style='white',font_scale=1.8)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
    figure = plt.figure(figsize=(3,3))
    ax = figure.add_subplot(111)
    ax.scatter(width_microns_eff, height_microns_eff, c='k', alpha=.45, s=100, lw=0)
    #ax.set_xlim(200,350)
    #ax.set_ylim(490, 340)
    xlim = ax.get_xlim()
    #ylim = ax.get_ylim()

    ax.set_xlim(.85*xlim[0], 1.15*xlim[1])
    #ax.set_ylim(.85*ylim[1], 1.15*ylim[0])

    ax.axvline(x=mean_width, ls=':')
    ax.axhline(y=mean_height,ls=':')

    ax.set_xlabel('Eff. FOV width [um]')
    ax.set_ylabel('Eff. FOV height [um]')
    figure.savefig(figure_dir / 'FOV across animals.pdf', bbox_inches='tight')


    return








if __name__ == "__main__":
   

    all_animals = [
                '90222','90218','90647',
                '82913','88592','89622',
                '87244','89841','60480',
                '87245','87187','88106',
                '94557','97045','97046',
                ] 

    animals_trans  = ['88592','82913','87187','88106','90222','90218','90647',
                      '89841','89622','94557','97045','97046']
    animals_virus  = ['87244','87245','60480']
    print(f'Including {len(all_animals)} animals')

    param_hash_id_cell ='standard'
    param_hash_session ='cf83e1357eefb8bd'

    # Fig. S1B
    print('\nFOV (eff) across animals')
    _get_fov_size(all_animals)
    
    # Fig. S1D
    print('\nSNR comparison transgenic vs. virus injected')
    _get_snr_virus_trans(animals_trans = animals_trans,
                         animals_virus = animals_virus,
                         param_hash_id_cell ='standard',
                         param_hash_session ='cf83e1357eefb8bd')

    # Fig. S1E
    print('\nCell count per session')
    _get_cell_count(all_animals, 
                    animals_trans = animals_trans,
                    animals_virus = animals_virus, 
                    param_hash_id_cell ='standard',
                    param_hash_session ='cf83e1357eefb8bd',
                    min_no_cells = 10)
                    

    print(figure_dir)
    print('Success.')