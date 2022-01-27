### Session level score cutoffs
 
import sys, os
import datajoint as dj
import numpy as np 

#### LOAD DATABASE #########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)
 

from .anatomical_distribution import FilteredCells, FilteredCellsParams

@imhotte
class CutoffsInfoContent(dj.Computed):    
    definition = """
    # Session level info content cutoffs 
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles               : int                  # Number of shuffles (total)
    info_content_90 = NULL   : double               # Spatial information content 90th cutoff
    info_content_95 = NULL   : double               # Spatial information content 95th cutoff
    info_content_99 = NULL   : double               # Spatial information content 99th cutoff 
    """
    
    @property
    def key_source(self):
        return FilteredSessions.proj() * FilteredCellsParams & Shuffled.RatemapStats 

    def make(self, key):
        # Session spatial information cutoffs
        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (Shuffled.RatemapStats & cell_keys & parameter_dict_cell).fetch('information_content_shuffles')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles']      = len(shuffles)
        key['info_content_90'] = np.nanpercentile(shuffles,90)
        key['info_content_95'] = np.nanpercentile(shuffles,95)
        key['info_content_99'] = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 

@imhotte
class CutoffsGridscore(dj.Computed):    
    definition = """
    # Session level gridscore cutoffs 
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles               : int                  # Number of shuffles (total)
    gridscore_90 = NULL      : double               # GridScore 90th cutoff
    gridscore_95 = NULL      : double               # GridScore 95th cutoff
    gridscore_99 = NULL      : double               # GridScore 99th cutoff 
    """

    @property
    def key_source(self):
        return FilteredSessions.proj() * FilteredCellsParams & Shuffled.GridScore 

    def make(self, key):

        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (Shuffled.GridScore & cell_keys & parameter_dict_cell).fetch('gridscore_shuffles')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles']   = len(shuffles)
        key['gridscore_90'] = np.nanpercentile(shuffles,90)
        key['gridscore_95'] = np.nanpercentile(shuffles,95)
        key['gridscore_99'] = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 


###### OVSCORES : SPECIAL CASE ! base_session instead of session_name 
@imhotte
class CutoffsOVScore(dj.Computed):    
    definition = """
    # Session level OV score cutoffs 
    -> FilteredSessions.proj(base_session='session_name')
    -> FilteredCellsParams
    ---
    n_shuffles               : int                  # Number of shuffles (total)
    ovscore_90 = NULL        : double               # OV Score 90th cutoff
    ovscore_95 = NULL        : double               # OV Score 95th cutoff
    ovscore_99 = NULL        : double               # OV Score 99th cutoff 
    """

    @property
    def key_source(self):
        return FilteredSessions.proj(base_session='session_name') * FilteredCellsParams & OVCScores # Only those that have OVCscores entries

    def make(self, key):

        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (OVCScores.ShuffledOVScore & cell_keys & parameter_dict_cell).fetch('shuffled_ovscores')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles'] = len(shuffles)
        key['ovscore_90'] = np.nanpercentile(shuffles,90)
        key['ovscore_95'] = np.nanpercentile(shuffles,95)
        key['ovscore_99'] = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 

####################################

@imhotte
class CutoffsMVL(dj.Computed):    
    definition = """
    # Session level MVL (head direction tuning) cutoffs 
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles               : int              # Number of shuffles (total)
    mvl_90 = NULL            : double           # MVL 90th cutoff
    mvl_95 = NULL            : double           # MVL 95th cutoff
    mvl_99 = NULL            : double           # MVL 99th cutoff 
    """

    @property
    def key_source(self):
        return FilteredSessions.proj() * FilteredCellsParams & Shuffled.AngularRateStats 

    def make(self, key):

        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (Shuffled.AngularRateStats & cell_keys & parameter_dict_cell).fetch('mvl_shuffles')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles'] = len(shuffles)
        key['mvl_90']     = np.nanpercentile(shuffles,90)
        key['mvl_95']     = np.nanpercentile(shuffles,95)
        key['mvl_99']     = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 


@imhotte
class CutoffsBorderScore(dj.Computed):    
    definition = """
    # Session level Borderscore cutoffs 
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles               : int              # Number of shuffles (total)
    borderscore_90 = NULL    : double           # Borderscore 90th cutoff
    borderscore_95 = NULL    : double           # Borderscore 95th cutoff
    borderscore_99 = NULL    : double           # Borderscore 99th cutoff 
    """

    @property
    def key_source(self):
        return FilteredSessions.proj() * FilteredCellsParams & Shuffled.BorderScore 

    def make(self, key):

        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (Shuffled.BorderScore & cell_keys & parameter_dict_cell).fetch('borderscore_shuffles')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles'] = len(shuffles)
        key['borderscore_90']     = np.nanpercentile(shuffles,90)
        key['borderscore_95']     = np.nanpercentile(shuffles,95)
        key['borderscore_99']     = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 


@imhotte
class CutoffsBVScore(dj.Computed):    
    definition = """
    # Session level boundary vector score cutoffs 
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_shuffles               : int              # Number of shuffles (total)
    bvs_90 = NULL            : double           # BVS 90th cutoff
    bvs_95 = NULL            : double           # BVS 95th cutoff
    bvs_99 = NULL            : double           # BVS 99th cutoff 
    """

    @property
    def key_source(self):
        return FilteredSessions.proj() * FilteredCellsParams & ShuffledBVS.BVS

    def make(self, key):

        cell_keys = FilteredCells & key
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        shuffles = (ShuffledBVS.BVS & cell_keys & parameter_dict_cell).fetch('bvs_shuffles')
        shuffles = np.concatenate(shuffles)

        key['n_shuffles'] = len(shuffles)
        key['bvs_90']     = np.nanpercentile(shuffles,90)
        key['bvs_95']     = np.nanpercentile(shuffles,95)
        key['bvs_99']     = np.nanpercentile(shuffles,99)

        self.insert1(key)
        return 




