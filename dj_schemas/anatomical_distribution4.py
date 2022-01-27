### Scores, etc. over anatomical space - second 
import numpy as np
import datajoint as dj
from pointpats import PointPattern, PoissonPointProcess
from sklearn.metrics import pairwise_distances

from tqdm.auto import tqdm

#### LOAD DATABASE #########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)
 

from .constants import GLOBAL_CUTOFF_NO_CELLS, GLOBAL_CUTOFF_NO_STARTER_CELLS, MIN_NN_DIST_THRESH

@imhotte 
class PairwDistParams(dj.Lookup):
    definition = """
    # Parameters for pairwise distance analysis
    pairwise_dist_param          : char(1)         # Param ID
    ---    
    score                        : varchar(50)     # Score (column) name
    score_cutoff                 : varchar(50)     # Score cutoff value (defines starter cell population)
    scoretables                  : varchar(1000)   # Datajoint tables that 'score' can be found in (e.g. 'GridScore')
    """


@imhotte 
class PairwDist(dj.Computed):
    definition = """
    # Pairwise distance statistic
    -> FilteredSessions
    -> FilteredCellsParams
    -> PairwDistParams
    -> ProjectionCorr
    ---    

    """

    class Cells(dj.Part):
        definition = """
        # Cell numbers
        -> master
        region                                  :  char(3)         # Brain region (3 letter abbreviation)
        --- 
        n_all                                   :  smallint        # How many cells in total were considered?
        n_startr                                :  smallint        # How many starter cells was the average calculated over (population A)? 
        """

    class PairwD(dj.Part):
        definition = """
        # Pairw. distance results per region
        -> master
        region                                  :  char(3)         # Brain region (3 letter abbreviation)
        ---
        med_pairw_dist                          :  double          # Median pairwise distance
        mean_pairw_dist                         :  double          # Mean pairwise distance
        med_pairw_dist_shuffall                 :  double          # Shuffled median pairwise distance taking all cells as start population
        mean_pairw_dist_shuffall                :  double          # Shuffled mean pairwise distance taking all cells as start population
        med_pairw_dist_shuffref = NULL          :  double          # Shuffled median pairwise distance taking only reference cells as start population
        mean_pairw_dist_shuffref = NULL         :  double          # Shuffled mean pairwise distance taking only reference cells as start population
        med_pairw_dist_csr                      :  double          # Shuffled median pairwise distance taking CSR as start population
        mean_pairw_dist_csr                     :  double          # Shuffled mean pairwise distance taking CSR as start population  
        """

    class NN(dj.Part):
        definition = """
        # Nearest neighbour (NN) distance results per region
        -> master
        region                                  :  char(3)         # Brain region (3 letter abbreviation)
        --- 
        mean_nn                                 :  blob@hottestore # Mean NN over 1-10 neighbours
        mean_nn_shuff_all                       :  blob@hottestore # Shuffled mean NN over 1-10 neighbours taking all cells as start population
        mean_nn_shuff_ref = NULL                :  blob@hottestore # Shuffled mean NN over 1-10 neighbours taking only reference cells as start population
        mean_nn_csr                             :  blob@hottestore # Shuffled mean NN over 1-10 neighbours taking CSR as start population  
        """

    @property
    def key_source(self):
        return super().key_source & RoisCorrBrainLoc # Only sessions that have RoisCorrBrainLoc() populated




    def make(self, key): 

        n_shuffles = 1000 # Fixed! 

        params = (PairwDistParams & key).fetch1()
        cell_parameter_dict = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        filtered_cells = eval(params['scoretables']) * RoisCorr & FilteredCells & key & cell_parameter_dict

        # GET MEC, PAS, ... ANATOMICAL FILTERS
        mec_filter = RoisCorrBrainLoc.MEC & key
        pas_filter = RoisCorrBrainLoc.PAS & key

        # Define entry dictionaries to loop over for part tables
        entry_dicts = [
                    {
                'cell_filter' : mec_filter,
                'name'        : 'MEC',
                    },
                    {
                'cell_filter' : pas_filter,
                'name'        : 'PAS',
                    },
                ]




        # Insert into master
        self.insert1(key, ignore_extra_fields=True)

        ### LOOP OVER BRAIN REGIONS 
        for entry_dict in entry_dicts:
            
            cells_session_part = filtered_cells & entry_dict['cell_filter'] # i.e. MEC cells, ...
            n_cells = len(cells_session_part)

            if (n_cells < GLOBAL_CUTOFF_NO_CELLS): 
                continue

            # Define starter cells
            starter_cells = cells_session_part & f'{params["score"]} > {params["score_cutoff"]}'
            n_starter_cells = len(starter_cells)
            if n_starter_cells < GLOBAL_CUTOFF_NO_STARTER_CELLS:
                continue

            cell_centers_all      = np.stack(cells_session_part.fetch('center_x_corr','center_y_corr')).T
            cell_centers_starters = np.stack(starter_cells.fetch('center_x_corr','center_y_corr')).T
            # Reference cells 
            ref_cell_ids = (cells_session_part - starter_cells.proj()).fetch('cell_id').tolist()
            cell_centers_ref      = np.stack((cells_session_part & [f'cell_id = {cell}' for cell in ref_cell_ids]).fetch('center_x_corr','center_y_corr')).T


            # DATA ...
            mean_nns, med_pairw_dist, mean_pairw_dist = get_nn_pwd(cell_centers=cell_centers_starters, 
                                                                   min_dist_thresh=MIN_NN_DIST_THRESH)


            # SHUFFLE ...
            # ALL CELLS 
            mean_nns_shuff_all, med_pairw_dist_shuff_all, mean_pairw_dist_shuff_all = get_shuffled_nn_pwd(
                                                                              cell_centers=cell_centers_all, 
                                                                              min_dist_thresh=MIN_NN_DIST_THRESH,  
                                                                              n_cells=n_starter_cells,
                                                                              n_shuffles=n_shuffles,
                                                                              points_csr='points'
                                                                              )
                

            # REFERENCE
            if len(cell_centers_ref) >= len(cell_centers_starters):
                mean_nns_shuff_ref, med_pairw_dist_shuff_ref, mean_pairw_dist_shuff_ref = get_shuffled_nn_pwd(
                                                                                cell_centers=cell_centers_ref, 
                                                                                min_dist_thresh=MIN_NN_DIST_THRESH,  
                                                                                n_cells=n_starter_cells,
                                                                                n_shuffles=n_shuffles,
                                                                                points_csr='points'
                                                                                )
            else: 
                mean_nns_shuff_ref, med_pairw_dist_shuff_ref, mean_pairw_dist_shuff_ref = None, None, None


            # CSR
            # Make a point pattern from starter cells and then create CSR over those 
            pp_starters = []
            for cell in zip(cell_centers_starters[:,0], cell_centers_starters[:,1]):
                pp_starters.append(cell)
            pp_starters = PointPattern(pp_starters)

            # Use PoissonPointProcess to create (uniform) CSR (complete spatial randomness)
            # https://pysal.org/notebooks/explore/pointpats/process.html
            csr_n = PoissonPointProcess(pp_starters.window, 
                                        pp_starters.n, 
                                        conditioning=False, 
                                        samples=n_shuffles, 
                                        asPP=True) 

            mean_nns_csr, med_pairw_dist_csr, mean_pairw_dist_csr = get_shuffled_nn_pwd(cell_centers=csr_n, 
                                                                                        min_dist_thresh=MIN_NN_DIST_THRESH,
                                                                                        points_csr='csr'
                                                                                        )


            # SAVE 
            entry_dict_celln = {
                 'region'                  : entry_dict['name'],
                 'n_all'                   : n_cells,
                 'n_startr'                : n_starter_cells,
            }

            entry_dict_pairwd = {
                 'region'                  : entry_dict['name'],
                 'med_pairw_dist'          : med_pairw_dist,
                 'mean_pairw_dist'         : mean_pairw_dist,
                 'med_pairw_dist_shuffall' : med_pairw_dist_shuff_all,
                 'mean_pairw_dist_shuffall': mean_pairw_dist_shuff_all,
                 'med_pairw_dist_shuffref' : med_pairw_dist_shuff_ref,
                 'mean_pairw_dist_shuffref': mean_pairw_dist_shuff_ref,
                 'med_pairw_dist_csr'      : med_pairw_dist_csr,
                 'mean_pairw_dist_csr'     : mean_pairw_dist_csr,                
            }

            entry_dict_nn = {
                 'region'                  : entry_dict['name'],
                 'mean_nn'                 : mean_nns,
                 'mean_nn_shuff_all'       : mean_nns_shuff_all,
                 'mean_nn_shuff_ref'       : mean_nns_shuff_ref,
                 'mean_nn_csr'             : mean_nns_csr,
            }


            self.Cells.insert1({**key, **entry_dict_celln}, ignore_extra_fields=True)
            self.PairwD.insert1({**key, **entry_dict_pairwd}, ignore_extra_fields=True)
            self.NN.insert1({**key, **entry_dict_nn}, ignore_extra_fields=True)








#### HELPERS FOR NEAREST NEIGHBOUR 

def get_nn_pwd(cell_centers, min_dist_thresh):
    '''
    Helper for pairwise distance analysis of 
    one (i.e. grid cell) population.
    
    Parameter
    ---------
    cell_centers   : np.array
                    2D array (n cells x dim) of cell centers
    
    min_dist_thresh : float
                    Minimum euclid. distance:
                    Only analyse distances above this threshold 
    
    Returns
    -------
    mean_nnd     : float
                Average nearest neighbour (NN) distance
    median_pwd   : float
                Median pairwise distances
    mean_pwd     : float
                Mean pairwise distances
                    
    '''
    
    # Get point pattern  
    pp = []
    for cell in zip(cell_centers[:,0], cell_centers[:,1]):
        pp.append(cell)
    pp = PointPattern(pp)
    
    # ... mean nearest neighbour distance
    # Over different number of nns 
    
    mean_nn_dict = {}
    for nn in range(1,11):
        _, nnd_ = pp.knn(nn)
        filtered_nnd_ = nnd_[nnd_ > min_dist_thresh]
        mean_nn_dict[nn] = np.nanmean(filtered_nnd_)
        
    # All pairwise distances 
    dist_mat = pairwise_distances(cell_centers, metric='euclidean')
    dist_mat = dist_mat[np.triu_indices_from(dist_mat, k=1)]
    filtered_dist_mat = dist_mat[dist_mat > min_dist_thresh]
    
    median_pairw_dist = np.nanmedian(filtered_dist_mat)
    mean_pairw_dist   = np.nanmean(filtered_dist_mat)
    
    return [mean_nn_dict[1], mean_nn_dict[2], mean_nn_dict[3], \
        mean_nn_dict[4], mean_nn_dict[5], mean_nn_dict[6], \
        mean_nn_dict[7], mean_nn_dict[8], mean_nn_dict[9], \
        mean_nn_dict[10]], \
        median_pairw_dist, mean_pairw_dist


def get_shuffled_nn_pwd(cell_centers, min_dist_thresh, n_cells=None, n_shuffles=1000,  points_csr='points'):
    '''
    Pick random cells from "cell_centers" of length "n_cells"
    OR 
    use existing PySAL CSR (complete spatial randomness) points.
    
    call get_nn_pwd() in the process to get NN / pairwise distance results.
    
    '''
    
    rnd_mean_nns1      = []
    rnd_mean_nns2      = []
    rnd_mean_nns3      = []
    rnd_mean_nns4      = []
    rnd_mean_nns5      = []
    rnd_mean_nns6      = []
    rnd_mean_nns7      = []
    rnd_mean_nns8      = []
    rnd_mean_nns9      = []
    rnd_mean_nns10     = []

    rnd_median_pairwd  = []
    rnd_mean_pairwd    = []
    
    
    
    if points_csr == 'points':
        if n_cells is None: 
            raise ValueError('Need to specify number of cells to pick')

        assert len(cell_centers) >= n_cells , f'# cells to pick ({n_cells}) is bigger than # all given cells'

        # Pick random points 
        for _ in tqdm(range(n_shuffles)):
            rnd_idxs  = np.random.choice(np.arange(len(cell_centers)), size=n_cells, replace=False)
            assert len(set(rnd_idxs)) == len(rnd_idxs)

            rnd_cells = cell_centers[rnd_idxs] 

            mean_nns, med_pairw_dist_, mean_pairw_dist_ = get_nn_pwd(rnd_cells, min_dist_thresh)

            rnd_mean_nns1.append(mean_nns[0])
            rnd_mean_nns2.append(mean_nns[1])
            rnd_mean_nns3.append(mean_nns[2])
            rnd_mean_nns4.append(mean_nns[3])
            rnd_mean_nns5.append(mean_nns[4])
            rnd_mean_nns6.append(mean_nns[5])
            rnd_mean_nns7.append(mean_nns[6])
            rnd_mean_nns8.append(mean_nns[7])
            rnd_mean_nns9.append(mean_nns[8])
            rnd_mean_nns10.append(mean_nns[9])

            rnd_median_pairwd.append(med_pairw_dist_)
            rnd_mean_pairwd.append(mean_pairw_dist_)
            
    elif points_csr=='csr': 
        if n_cells is not None:
            raise ValueError('n_cells is not None: check parameters')
        # Loop over realizations of CSR 
        for i in tqdm(range(cell_centers.samples), desc='CSR'):
            points_ = cell_centers.realizations[i].points.values
            
            mean_nns, med_pairw_dist_, mean_pairw_dist_ = get_nn_pwd(points_, min_dist_thresh)

            rnd_mean_nns1.append(mean_nns[0])
            rnd_mean_nns2.append(mean_nns[1])
            rnd_mean_nns3.append(mean_nns[2])
            rnd_mean_nns4.append(mean_nns[3])
            rnd_mean_nns5.append(mean_nns[4])
            rnd_mean_nns6.append(mean_nns[5])
            rnd_mean_nns7.append(mean_nns[6])
            rnd_mean_nns8.append(mean_nns[7])
            rnd_mean_nns9.append(mean_nns[8])
            rnd_mean_nns10.append(mean_nns[9])

            rnd_median_pairwd.append(med_pairw_dist_)
            rnd_mean_pairwd.append(mean_pairw_dist_)
        
    else:
        raise NotImplementedError(f'{points_csr} not implemented')
        
    return [np.nanmedian(rnd_mean_nns1), np.nanmedian(rnd_mean_nns2), np.nanmedian(rnd_mean_nns3),\
           np.nanmedian(rnd_mean_nns4), np.nanmedian(rnd_mean_nns5), np.nanmedian(rnd_mean_nns6),\
           np.nanmedian(rnd_mean_nns7), np.nanmedian(rnd_mean_nns8), np.nanmedian(rnd_mean_nns9),\
           np.nanmedian(rnd_mean_nns10)], \
           np.nanmedian(rnd_median_pairwd), np.nanmedian(rnd_mean_pairwd)