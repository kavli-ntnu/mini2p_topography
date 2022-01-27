### Helpers for pairwise distance notebook
from pathlib import Path
import numpy as np
from scipy.stats import wasserstein_distance
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns 

# To prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

#### NEW FUNCTIONS MARCH 2021

def norm_pairw_nn_df(df, cols_to_norm, cols, norm_to):
    '''
    Return normalized pairwise distance or NN results
    (Fetched from PairwDist.PairwD or PairwDist.NN)

    Parameter
    ---------
    df : Pandas dataframe with PairwDist results
    cols_to_norm: list: Columns that should be kept and normalized 
    cols : list : Additional columns that should be kept but not normalized
    norm_to : string: Column name that the others (cols) should be normalized to
    '''
    if isinstance(cols_to_norm, str):
        cols_to_norm = [cols_to_norm]
    if isinstance(cols, str):
        cols = [cols]
    if isinstance(norm_to, str):
        norm_to = [norm_to]
    
    print(f'Normalising to {norm_to}')
    
    df_kept = df.copy()
    df_kept = df_kept[set(cols_to_norm + cols + norm_to)]

    # Normalize
    for col in df_kept.columns.values:       
        if col in norm_to + cols:
            continue
        else:
            df_kept[col] = df_kept[col] / df_kept[norm_to[0]]
    df_kept[norm_to] = np.ones_like(df_kept[norm_to])
    
    return df_kept


def norm_nn_df(df, norm_to='mean_nn_shuff_all'):
    '''
    Return normalized NN result (NN over n nearest neighbours)

    Parameter
    ---------
    df : Pandas dataframe with PairwDist.NN() results
    norm_to : string: Column name that "mean_nn" (data) should be normalized to
                      CAVE Only possible option atm is norm_to = mean_nn_shuff_all

    Returns
    -------
    mean_nns :            numpy array
                          Normalized mean NN distances over NN 
    mean_nn_shuff_refs :  numpy array
                          Normalized mean NN distances of reference population over NN 
    '''
    assert norm_to == 'mean_nn_shuff_all', f'Normalisation of NN to {norm_to} is not implemented yet'
    print(f'Normalising mean_nn and mean_nn_shuff_ref to {norm_to}')
    
    mean_nns = []
    mean_nn_shuff_refs = []

    # Normalize all to 'mean_nn_shuff_all'
    for _, nn in df.iterrows():
        mean_nns.append(np.array(nn['mean_nn']) / np.array(nn[norm_to]))
        mean_nn_shuff_refs.append(np.array(nn['mean_nn_shuff_ref']) / np.array(nn[norm_to]))

    mean_nns = np.stack(mean_nns)
    mean_nn_shuff_refs = np.stack(mean_nn_shuff_refs)
    return mean_nns, mean_nn_shuff_refs






def plot_pairw_nn_summary(pairw_df_norm, 
                          cols_to_norm, 
                          colors, 
                          xlabels=None, 
                          save_path=None, 
                          label=''
                          ):
    '''
    Plot pairwise distance or NN summary (line + boxplots)
    '''

    sns.set(style='white',font_scale=1.4)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True

    figure = plt.figure(constrained_layout=False, figsize=(3,4))
    gs = figure.add_gridspec(nrows=1, ncols=9, wspace=0, hspace=1)
    ax_lines   =  figure.add_subplot(gs[:, :5])
    ax_box     =  figure.add_subplot(gs[:, 6:]) 

    no_entries      = len(pairw_df_norm)
    sqrt_no_entries = np.sqrt(no_entries)

    # LINE PLOT
    for no, row in pairw_df_norm.iterrows():
        ax_lines.plot(np.arange(len(cols_to_norm)),
                np.array([row[col] for col in cols_to_norm]),
                color=colors[no], alpha=.3)

    averages = np.array([np.nanmean(pairw_df_norm[col]) for col in cols_to_norm])
    sems     = np.array([np.std(pairw_df_norm[col])/sqrt_no_entries for col in cols_to_norm])

    ax_lines.errorbar(np.arange(len(cols_to_norm)),
                      averages,
                      sems,
                      color='k', lw=3, zorder=10, marker='.',alpha=.9
                      )

    ax_lines.set_xticks(np.arange(len(cols_to_norm)))
    if xlabels is None:
        ax_lines.set_xticklabels([col for col in cols_to_norm], rotation=40,ha='right')
    else:
        ax_lines.set_xticklabels(xlabels, rotation=40, ha='right')
    ax_lines.axhline(y=1, ls=':', color='k')
    
    ax_lines.set_xlim([-.25, len(cols_to_norm)-.75])
    ax_lines.set_ylabel(f'Norm. {label} distance')
    ylim = ax_lines.get_ylim()
    
    #BOX PLOT
    ax_box.boxplot(pairw_df_norm[cols_to_norm], widths=.7, showmeans=True, meanprops={'marker':'+','markeredgecolor':'k'})
    ax_box.axhline(y=1, ls=':', color='k')
    ax_box.set_xticks(np.arange(len(cols_to_norm))+1)
    if xlabels is None:
        ax_box.set_xticklabels([col for col in cols_to_norm], rotation=40, ha='right')
    else:
        ax_box.set_xticklabels(xlabels, rotation=40, ha='right')
    ax_box.set_ylim(ylim[0], ylim[1])
    ax_box.get_yaxis().set_visible(False)
    
    sns.despine(left=True)

    if save_path is not None: 
        save_path = Path(save_path)
        figure.savefig(save_path / f'dist_summary_{label}.pdf', dpi=300, bbox_inches='tight')

    #plt.show()



def plot_mean_nn_over_nn(mean_nns, mean_nn_shuff_refs, save_path=None):
    '''
    PairwDist.NN()
    Create plot of normalised mean NN distance over number
    of nearest neighbours (PairwDist.NN)
    
    
    '''
    sns.set(style='white',font_scale=1.5)
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True



    figure = plt.figure(figsize=(3,4))
    ax = figure.add_subplot(111)

    ax.axhline(y=1, ls=':', color='k', lw=1)

    assert mean_nns.shape == mean_nn_shuff_refs.shape, 'Matrices do not have the same dimensions'

    # Plot mean / std
    n = mean_nns.shape[0]
    sqrt_n = np.sqrt(n)

    # Ref
    for row in np.arange(mean_nn_shuff_refs.shape[0]):
        ax.plot(np.arange(mean_nn_shuff_refs.shape[1]),
                mean_nn_shuff_refs[row,:],
                color='cornflowerblue', lw=1, alpha=.1)


    ax.errorbar(np.arange(mean_nn_shuff_refs.shape[1]),
                    np.nanmean(mean_nn_shuff_refs, axis=0),
                    np.nanstd(mean_nn_shuff_refs, axis=0)/sqrt_n,
                    color='cornflowerblue', lw=2, zorder=10, marker='.',alpha=.9
                    )

    # Data
    for row in np.arange(mean_nns.shape[0]):
        ax.plot(np.arange(mean_nns.shape[1]),
                mean_nns[row,:],
                color='k', lw=1, alpha=.1)


    ax.errorbar(np.arange(mean_nns.shape[1]),
                    np.nanmean(mean_nns, axis=0),
                    np.nanstd(mean_nns, axis=0)/sqrt_n,
                    color='k', lw=2, zorder=10, marker='.',alpha=.9
                    )

    ax.set_xticks(np.arange(mean_nns.shape[1]))
    ax.set_xticklabels(np.arange(mean_nns.shape[1])+1)
    ax.set_xlim([-.25, mean_nns.shape[1]-.75])

    sns.despine(left=True,bottom=True)

    if save_path is not None: 
        save_path = Path(save_path)
        figure.savefig(save_path / f'mean_nn_over_nn.pdf', dpi=300, bbox_inches='tight')



