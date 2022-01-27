### Helpers for pairwise distance notebook
from pathlib import Path
import numpy as np
from tqdm.auto import tqdm
from sklearn.neighbors import KDTree
from pointpats import PointPattern, PoissonPointProcess

# Plotting
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white')

# To prevent bug in figure export as pdf: 
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


#### SIMULATIONS / TESTS ########################################################################################
# DATA (with subsampling)

def _mean_pp_nn(pp, min_dist_thresh):
    # Given a pysal pointpattern (pp), return mean nearest neighbour distance above 
    # minium distance threshold (min_dist_thresh)
    # https://github.com/pysal/pointpats/blob/0ca9085fdb082194540568166343b97111993272/pointpats/pointpattern.py#L324
    _, nnd = pp.knn(1)
    filtered_nnd = nnd[nnd > min_dist_thresh]
    mean_nnd = np.nanmean(filtered_nnd)
    return mean_nnd



def calc_nn(pop_A, pop_B, min_dist_thresh):
    '''
    Given center positions in two populations (pop_A vs. pop_B),
    return nearest neighbour distances (NN) A->B and B->A.

    Calculate inter to intra population NN distances ratio.

    '''

    # ... build tree
    # leafsize: The number of points at which the algorithm switches over to brute-force. Default: 10
    tree_A = KDTree(pop_A, leaf_size=10, metric='euclidean') 
    tree_B = KDTree(pop_B, leaf_size=10, metric='euclidean')
    
    NN_AB = [] # Nearest neighbour distances A->B
    NN_BA = [] # ......                      B->A
    
    # ... loop over points
    points_ = []
    for point in pop_A: 
        points_.append(point)

        # For each point in population A ... 
        dist, _ = tree_B.query([point], k=1)
        dist = dist[0][0]
        if dist > min_dist_thresh:
            NN_AB.append(dist)
    points_ = PointPattern(points_)
    # Ratio inter to intra A->B
    ratio_AB = (np.nanmean(NN_AB) / _mean_pp_nn(points_, min_dist_thresh))

    points_ = []
    for point in pop_B: 
        points_.append(point)

        # For each point in population B ... 
        dist, _ = tree_A.query([point], k=1)   
        dist = dist[0][0]
        if dist > min_dist_thresh:
            NN_BA.append(dist)
    points_ = PointPattern(points_) 
    # Ratio inter to intra B->A
    ratio_BA = (np.nanmean(NN_BA) / _mean_pp_nn(points_, min_dist_thresh))

    return NN_AB, NN_BA, ratio_AB, ratio_BA




def get_nn_pwd_dists(pop_A, pop_B, min_dist_thresh, chunk_length=None, permutations=1000): 
    '''
    
    For two populations A and B, get
    - NN distances for A vs. B and B vs. A 
    - Ratio of inter vs. intra cluster distances

    Parameters
    ----------
    pop_A   : np.array 
              Array of centers (number of cells x dim)
    pop_B   : np.array 
              Array of centers (number of cells x dim)
    min_dist_thresh : float
                      Minimum distance threshold
    chunk_length : None or int
                   If None, use full size, if int: Cut out sub-population
                   of chunk_length 
    permutations : int
                   Number of permutations
    '''
    
 

    if chunk_length is None: 
        # Take "full" population size (do not sub-divide)
        NN_AB, NN_BA, ratio_AB, ratio_BA   = calc_nn(pop_A, pop_B, min_dist_thresh)
        pop_A_sub_sample, pop_B_sub_sample = pop_A, pop_B

    else: 
        # Permute, cut out chunk length 
        nnsAB = [] 
        nnsBA = []        
        ratios_AB = [] 
        ratios_BA = []
        
        for iter in tqdm(range(permutations), desc='NN data'):
            pop_A = np.random.permutation(pop_A)
            pop_B = np.random.permutation(pop_B)

            pop_A_sub = pop_A[:chunk_length]
            pop_B_sub = pop_B[:chunk_length]
            assert len(pop_A_sub) == chunk_length, f'Population size A differs from chunk length ({len(pop_A_sub)} | {chunk_length})'
            assert len(pop_B_sub) == chunk_length, f'Population size B differs from chunk length ({len(pop_B_sub)} | {chunk_length})'
            # Save example 
            if not iter: 
                pop_A_sub_sample, pop_B_sub_sample = pop_A_sub, pop_B_sub 

            NN_AB, NN_BA, ratio_AB, ratio_BA = calc_nn(pop_A_sub, pop_B_sub, min_dist_thresh)
            nnsAB.append(NN_AB)
            nnsBA.append(NN_BA)
            ratios_AB.append(ratio_AB)
            ratios_BA.append(ratio_BA)

        # plt.hist(ratios_AB,lw=0)
        # plt.show() 
        # plt.hist(ratios_BA,lw=0)
        # plt.show()

        NN_AB       = np.array([item for sublist in nnsAB for item in sublist])
        NN_BA       = np.array([item for sublist in nnsBA for item in sublist])
        ratio_AB    = np.nanmedian(ratios_AB)
        ratio_BA    = np.nanmedian(ratios_BA)

    # Save output 
    nns = {
        'AB' : NN_AB,
        'BA' : NN_BA
    }
    ratios = {
        'AB' : ratio_AB,
        'BA' : ratio_BA,
    }
    
    return nns, ratios, \
           pop_A_sub_sample, pop_B_sub_sample




def shuffle_nn_pwd_dists(centers, 
                         min_dist_thresh,
                         chunk_length=None, 
                         length_A=None, 
                         length_B=None, 
                         shuffle_iter=1000):
    '''
    
    Given a list of centers and a length of points for each chunk, 
    or the population size of A and B, 
    draw samples of len(chunk) or len(A) and len(B) 
    to generate two populations A and B, and 
    calculate NN statistics as in 
    `get_nn_pwd_dists()`.

    '''

    nnsAB_shuff = [] 
    nnsBA_shuff = [] 
       
    ratios_AB_shuff = [] 
    ratios_BA_shuff = []


    if (length_A is not None) and (length_B is not None):
        #print('Using population length instead of chunk length')
        use_pop_lengths = True
    else:
        use_pop_lengths = False

    for _ in range(shuffle_iter):
        
        if use_pop_lengths:
            # Get random indices of length len(population A) + len(population B)
            rnd_idxs_AB = np.random.choice(np.arange(len(centers)), size=length_A+length_B, replace=False)
            pop_A_      = centers[rnd_idxs_AB[:length_A]] # ... the rest is B
            pop_B_      = centers[rnd_idxs_AB[length_A:]]
        else: 
            # Get random indices twice the length of "chunk_length" and split in half
            # for the two "random populations"
            rnd_idxs    = np.random.choice(np.arange(len(centers)), size=chunk_length*2, replace=False)
            pop_A_      = centers[rnd_idxs[:chunk_length]]
            pop_B_      = centers[rnd_idxs[chunk_length:]]
        
        # # 1. Pairwise distances pop A vs. pop B    
        # SKIP FOR NOW
        # pop_dist = pairwise_distances(pop_A_,pop_B_).flatten()
        # pop_dists_shuff.append(pop_dist)

        # 2. NN         
        NN_AB, NN_BA, ratio_AB, ratio_BA = calc_nn(pop_A_, pop_B_, min_dist_thresh)

        nnsAB_shuff.append(NN_AB)
        nnsBA_shuff.append(NN_BA)
        ratios_AB_shuff.append(ratio_AB)
        ratios_BA_shuff.append(ratio_BA)

    
    # RESULTS
    # Raw NN distribution     
    nnsAB_shuff = np.array([item for sublist in nnsAB_shuff for item in sublist])
    nnsBA_shuff = np.array([item for sublist in nnsBA_shuff for item in sublist])
    
    nns_shuff = {
            'AB' : nnsAB_shuff,
            'BA' : nnsBA_shuff
        }

    # Ratios 
    median_ratio_inter_intra_AB_shuff = np.nanmedian(ratios_AB_shuff)
    median_ratio_inter_intra_BA_shuff = np.nanmedian(ratios_BA_shuff)

    ratios_shuff = {
            'AB' : median_ratio_inter_intra_AB_shuff,
            'BA' : median_ratio_inter_intra_BA_shuff,
        }
    
    return nns_shuff, ratios_shuff



def csr_nn_pwd_dists(cell_centers_A, 
                     cell_centers_B, 
                     min_dist_thresh,
                     chunk_length=None, 
                     length_A=None, 
                     length_B=None,  
                     shuffle_iter=1000):
    '''
    Given two populations A and B and their centers, 
    generate a CSR (complete spatial randomness) pattern. 


    '''
    # # Point process CSR over cell positions (window) in population A and B
    points = PointPattern(np.concatenate((cell_centers_A,cell_centers_B))) 
    csr_n  = PoissonPointProcess(points.window, points.n, conditioning=False, samples=shuffle_iter, asPP=True) 

    nn_csr_ab, nn_csr_ba = [], []
    ratios_csr_ab, ratios_csr_ba = [], [] 

    for i in range(shuffle_iter): # Create that many distinct CSR patterns
        points_ = csr_n.realizations[i].points.values    
        nns_csr, ratios_csr = shuffle_nn_pwd_dists(points_, 
                                                   min_dist_thresh, 
                                                   length_A=length_A, 
                                                   length_B=length_B, 
                                                   chunk_length=chunk_length,
                                                   shuffle_iter=1)
        nn_csr_ab.append(nns_csr['AB'])
        nn_csr_ba.append(nns_csr['BA'])
        ratios_csr_ab.append(ratios_csr['AB'])
        ratios_csr_ba.append(ratios_csr['BA'])       

    nn_csr_ab = np.array([item for sublist in nn_csr_ab for item in sublist])
    nn_csr_ba = np.array([item for sublist in nn_csr_ba for item in sublist])

    nns_csr = {
        'AB' : nn_csr_ab, 
        'BA' : nn_csr_ba,
    }

    ratios_csr = {
        'AB' : np.nanmedian(ratios_csr_ab),
        'BA' : np.nanmedian(ratios_csr_ba),
    }

    return nns_csr, ratios_csr




##### PLOTTING HELPERS #####################################################################################

def plot_blob_set(centers, 
                  labels, 
                  centroids=None, 
                  pop_A_sub=None, 
                  pop_B_sub=None, 
                  title=None, 
                  legend=True,
                  label_A='A',
                  label_B='B',
                  scalebar=50.):
    # Scalebar should be in microns if ProjectionCorr entries are used

    #colors = ['#979797', '#ee8c86']
    colors = ['#979797', '#75A2EA']

    if (pop_A_sub is not None) and (pop_B_sub is not None):
        sub = True
        s = 200
    else:
        sub = False
        s = 100

    sns.set(style='white',font_scale=1.3)
    figure = plt.figure(figsize=(6,6))
    ax = figure.add_subplot(111)

    for cluster in np.unique(labels):
        ax.scatter(centers[labels==cluster][:,0], 
                   centers[labels==cluster][:,1], 
                   color=colors[cluster], 
                   s=s, 
                   alpha=[.3 if sub else .85][0], 
                   lw=0)
        
    if sub:
        # Draw subsampled populations: 
        ax.scatter(pop_A_sub[:,0], 
                pop_A_sub[:,1], 
                color=colors[0], 
                s=40, 
                alpha=1, 
                lw=2,
                label=label_A)
        ax.scatter(pop_B_sub[:,0], 
                pop_B_sub[:,1], 
                color=colors[1], 
                s=40, 
                alpha=1, 
                lw=2,
                label=label_B)
        if legend:
            ax.legend()

    if centroids is not None: 
        for no, center in enumerate(centroids): 
            ax.scatter(center[0], center[1], marker='o', c='w', s=300, alpha=.8)
            ax.scatter(center[0], center[1], marker='x', c=[colors[0] if no==0 else colors[1]][0], s=300)


    XLIM = ax.get_xlim()
    YLIM = ax.get_ylim() 

    if title is not None:
        ax.set_title(title)
    
    plt.axis('scaled')

    # # Expand plotting area a bit (10%)
    # xlims = ax.get_xlim()
    # dxlims = np.abs(np.diff(xlims))
    # perc_dxlims = .1* dxlims
    # ylims = ax.get_ylim()
    # dylims = np.abs(np.diff(ylims))
    # perc_dylims = .1 * dylims
    # ax.set_xlim(xlims[0]-perc_dxlims, xlims[1]+perc_dxlims)
    # ax.set_ylim(ylims[0]-perc_dylims, ylims[1]+perc_dylims)

    ax.get_xaxis().set_ticks([]); ax.get_yaxis().set_ticks([]); 
    ax.invert_yaxis()

    if scalebar > 0:
         ax.plot([np.max(XLIM)-scalebar-2,
                  np.max(XLIM)-scalebar-2+scalebar], 
                 [np.max(YLIM)-2,
                  np.max(YLIM)-2], 
                  lw=4, color='k', 
                  alpha=.8, solid_capstyle='butt')


    sns.despine(left=True,bottom=True)
    #plt.show()

    return figure