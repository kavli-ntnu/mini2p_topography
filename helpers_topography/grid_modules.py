# Module finding
import numpy as np
import pandas as pd

from sklearn import preprocessing, neighbors
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, estimate_bandwidth
from hdbscan import HDBSCAN

# Plotting 
from dj_plotter.helpers.plotting_helpers import make_linear_colormap
# Make plots pretty 
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')


def extract_grid_modules(grid_stats, columns, 
                         scaling='MaxAbsScaler', 
                         reducing='pca', 
                         clustering='hdbscan', 
                         plot=True, 
                         save_path=None
                         ):
    '''
    Given a dataframe of grid statistics, find clusters (aka modules)
    
    Recommended stats:
    ['gs_orientation', 'gs_spacing', 'gs_ellipse_theta', 'gs_ellipse_aspect_ratio', 'avg_fieldsize']
    
    '''
    # Reduce dataframe to only essential (cluster worthy) columns
    grid_stats = grid_stats[columns]    
    print('Running stats over ... ')
    for column in columns:
        print(column)
    print('\n')
    #######################################################################################################################
    # Scaling / Normalization 
    if scaling == 'MaxAbsScaler': 
        print('Using MaxAbsScaler')
        scaler = preprocessing.MaxAbsScaler().fit(grid_stats)
        grid_stats_norm = scaler.transform(grid_stats)
    elif scaling == 'StandardScaler':
        print('Using StandardScaler')
        scaler = preprocessing.StandardScaler().fit(grid_stats)
        grid_stats_norm = scaler.transform(grid_stats)
    elif scaling is None:
        print('Skipping scaling')
        scaler = None
        grid_stats_norm = grid_stats
    else:
        raise NotImplementedError(f'{scaling} not implemented yet')


    #######################################################################################################################
    # Reduction
    if reducing == 'pca':
        print('Using PCA')
        reducer = PCA(n_components=2, whiten=True)
        reducer.fit(grid_stats_norm)
        reduced_gridstats = reducer.transform(grid_stats_norm)
        print(f'Explained variance by 2 components: {np.sum(reducer.explained_variance_ratio_):.2f}')
    elif reducing is None:
        print('Skipping decomposition')
        reducer = None
        reduced_gridstats = grid_stats_norm
    else:
        raise NotImplementedError(f'{reducing} not implemented yet')

    print(f'Data has now the following dimensions: {reduced_gridstats.shape}')


    #######################################################################################################################
    # Clustering    
    if clustering == 'meanshift':
        print('Using meanshift clustering')
        # Estimate the bandwidth to use with the mean-shift algorithm, and cluster
        bandwidth = estimate_bandwidth(reduced_gridstats, quantile=.2)
        clusterer = MeanShift(bandwidth=bandwidth, bin_seeding=True, min_bin_freq=2, cluster_all=False, max_iter=10000)
        clusterer.fit(reduced_gridstats)
    elif clustering == 'hdbscan':
        print('Using HDBSCAN clustering')
        clusterer = HDBSCAN(min_cluster_size=3, min_samples=2, cluster_selection_epsilon=.8, 
                            metric='manhattan', allow_single_cluster=False, prediction_data=True)
        clusterer.fit(reduced_gridstats)
    else:
        raise NotImplementedError(f'{clustering} not implemented yet')


    #######################################################################################################################
    # Labels and centers 
    # ... labels
    labels = clusterer.labels_
    # ... quick info
    labels_unique = np.unique(labels)
    n_clusters    = len(labels_unique[labels_unique>-1])
    print(f'Found {n_clusters} modules ({labels_unique})')
    
    # ... cluster centers
    if clustering != 'hdbscan':
        cluster_centers = clusterer.cluster_centers_
    else: # HDBSCAN needs special treatment
        cluster_centers = []
        for cluster in labels_unique[labels_unique > -1]:
            cluster_centers.append(clusterer.weighted_cluster_centroid(cluster))

    #######################################################################################################################
    # Classifier 
    classifier = neighbors.KNeighborsClassifier(5, metric='minkowski', algorithm='auto', weights='uniform')
    classifier.fit(reduced_gridstats, labels)

    #######################################################################################################################
    # Create figure? 

    if plot:
        colors = make_linear_colormap(np.arange(n_clusters), cmap='cmr.guppy')
        sns.set(style='white', font_scale=1.4)
        figure=plt.figure(figsize=(9,4))

        ax=figure.add_subplot(121)
        ax.scatter(reduced_gridstats[:,0], reduced_gridstats[:,1], s=5, color='k', alpha=.8)
        for no,center in enumerate(cluster_centers):
            ax.scatter(reduced_gridstats[:,0][labels==no], reduced_gridstats[:,1][labels==no], 
                        s=55, lw=0, color=colors[no], alpha=.7)
            ax.scatter(center[0],center[1], s=100, marker='x', color='k')
        ax.set_xlabel('C1')
        ax.set_ylabel('C2')
        ax.set_title(f'{reducing} + clustering')
        
        # Transfer labels to original spacing / orientation plot and draw 
        ax=figure.add_subplot(122)
        saved_figs_ = []
        ax.scatter(grid_stats.gs_orientation, grid_stats.gs_spacing, s=5, color='k', alpha=.8)
        for no,center in enumerate(cluster_centers):
            p = ax.scatter(grid_stats.gs_orientation.iloc[labels==no], grid_stats.gs_spacing.iloc[labels==no], 
                            color=colors[no], s=45, lw=0, alpha=.7, label=no)
            saved_figs_.append(p)
        ax.legend(handles=saved_figs_, title='Module', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.set_xlabel('Orientation')
        ax.set_ylabel('Spacing')
        ax.set_xlim(-35,35)

        plt.tight_layout()
        if save_path is not None: 
            figure.savefig(save_path + 'pca and clustering.pdf', bbox_inches='tight')
        plt.show()

    
    return scaler, reducer, clusterer, cluster_centers, labels, reduced_gridstats, classifier


def embed_existing_module(grid_stats, scaler, reducer, classifier):
    '''
    For a single (new) gridstats entry, use previously fitted 
    scaler, reducer and clusterer to find cluster label 

    Parameter
    ---------
    grid_stats     : pd.DataFrame or pd.Series 
    scaler         : (Fitted) sklearn.preprocessing instance
    reducer        : (Fitted) sklearn.decomposition instance
    classifier     : (Fitted) classifier 

    Returns
    -------
    reduced_gridstats : np.array
                        Scaled and reduced grid stats
    labels            : np.array
                        Class labels 

    '''

    # Test format
    if isinstance(grid_stats, pd.Series): # i.e. a single entry
        grid_stats = grid_stats.values.reshape(1, -1)

    if scaler is not None: 
        grid_stats_norm = scaler.transform(grid_stats)
    else:
        grid_stats_norm = grid_stats 

    if reducer is not None: 
        reduced_gridstats = reducer.transform(grid_stats_norm)
    else:
        reduced_gridstats = grid_stats_norm

    labels = classifier.predict(reduced_gridstats)

    return reduced_gridstats, labels.astype(int)


#### Plotting helper 
def plot_maps(maps, title='0', titles=None):
    '''
    Plot array of square maps
    
    '''
    sns.set(style='white', font_scale=1.2)

    if len(maps) > 81:
        print(f'Warning! Plotting at max 81 maps (given: {len(maps)})')
    figure = plt.figure(figsize=(10,10))
    for no,map_ in enumerate(maps): 
        ax=figure.add_subplot(9,9,no+1)
        ax.imshow(map_, cmap='cmr.neutral')
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        if titles is not None:
            ax.text(10,10, f'{titles[no]}', color='w')

        sns.despine(left=True,bottom=True)
        if no==80: break
    plt.suptitle(f'{title}', y=.92)
    #plt.tight_layout()
    plt.show()