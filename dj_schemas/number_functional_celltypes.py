### Number of functional cell types across sessions
 
from re import S
import sys, os
import datajoint as dj
import numpy as np 

#### LOAD DATABASE #########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)
 
from .anatomical_distribution import FilteredCells, FilteredCellsParams

@imhotte
class NumberCellTypes(dj.Computed):    
    definition = """
    # Number of (pure) cell types per session and brain region (chose MEC and PAS)
    -> FilteredSessions
    -> FilteredCellsParams
    ---
    n_cells                     : int          # Total number of (filtered) cells in this session
    """

    class MEC(dj.Part):
        definition = """
        # Cell numbers MEC
        -> master
        --- 
        mec_n_cells              : int         # MEC Total number of (filtered) cells
        mec_ovc = NULL           : int         # MEC Number of OV cells
        mec_n_grid_95            : int         # MEC Number of grid cells > 95th
        mec_n_info_95            : int         # MEC Number of cells w info content > 95th
        mec_n_hd_95              : int         # MEC Number of HD cells > 95th
        mec_n_border_95          : int         # MEC Number of Border cells > 95th
        mec_n_bv_95              : int         # MEC Number of BV cells > 95th
        mec_n_grid_99            : int         # MEC Number of grid cells > 99th
        mec_n_info_99            : int         # MEC Number of cells w info content > 99th
        mec_n_hd_99              : int         # MEC Number of HD cells > 99th
        mec_n_border_99          : int         # MEC Number of Border cells > 99th
        mec_n_bv_99              : int         # MEC Number of BV cells > 99th
        """

    class PAS(dj.Part):
        definition = """
        # Cell numbers PAS
        -> master
        --- 
        pas_n_cells              : int          # PAS Total number of (filtered) cells
        pas_ovc = NULL           : int          # PAS Number of OV cells
        pas_n_grid_95            : int          # PAS Number of grid cells > 95th
        pas_n_info_95            : int          # PAS Number of cells w info content > 95th
        pas_n_hd_95              : int          # PAS Number of HD cells > 95th
        pas_n_border_95          : int          # PAS Number of Border cells > 95th
        pas_n_bv_95              : int          # PAS Number of BV cells > 95th
        pas_n_grid_99            : int          # PAS Number of grid cells > 99th
        pas_n_info_99            : int          # PAS Number of cells w info content > 99th
        pas_n_hd_99              : int          # PAS Number of HD cells > 99th
        pas_n_border_99          : int          # PAS Number of Border cells > 99th
        pas_n_bv_99              : int          # PAS Number of BV cells > 99th
        """

    @property
    def key_source(self):
        return super().key_source & RoisCorrBrainLoc # Filter out those that do not have annotated brain regions


    def make(self, key):
        # Retrieve number of functional cell types per session
        cell_keys_session = FilteredCells & key
        # Brain region filter 
        cells_session_mec = cell_keys_session & RoisCorrBrainLoc.MEC
        cells_session_pas = cell_keys_session & RoisCorrBrainLoc.PAS

        # Cell parameter dictionary
        parameter_dict_cell = (FilteredCellsParams & key).fetch1('parameter_dict_cell')

        # Filter cells by score / cutoff 
        ovcells       = OVC.proj(is_ovc='is_ovc', session_name='base_session') & 'is_ovc > 0.5'

        bordercells95 = BorderScore * CutoffsBorderScore.proj('borderscore_95') & 'borderscore > borderscore_95'
        bvcells95     = BVScore * CutoffsBVScore.proj('bvs_95') & 'bvs > bvs_95'
        gridcells95   = GridScore * CutoffsGridscore.proj('gridscore_95') & 'gridscore > gridscore_95'
        hdcells95     = AngularRate.Stats * CutoffsMVL.proj('mvl_95') & 'mvl > mvl_95'
        infocells95   = Ratemap.Stats * CutoffsInfoContent.proj('info_content_95') & 'information_content > info_content_95'

        bordercells99 = BorderScore * CutoffsBorderScore.proj('borderscore_99') & 'borderscore > borderscore_99'
        bvcells99     = BVScore * CutoffsBVScore.proj('bvs_99') & 'bvs > bvs_99'
        gridcells99   = GridScore * CutoffsGridscore.proj('gridscore_99') & 'gridscore > gridscore_99'
        hdcells99     = AngularRate.Stats * CutoffsMVL.proj('mvl_99') & 'mvl > mvl_99'
        infocells99   = Ratemap.Stats * CutoffsInfoContent.proj('info_content_99') & 'information_content > info_content_99'


        ##### MEC ######################################################################################################################
        entry_dict_mec = {
            **key,
            'mec_n_cells'      : len(cells_session_mec & parameter_dict_cell),
            'mec_ovc'          : len(ovcells & cells_session_mec & parameter_dict_cell),
            'mec_n_grid_95'    : len(gridcells95 & cells_session_mec & parameter_dict_cell),
            'mec_n_info_95'    : len(infocells95 & cells_session_mec & parameter_dict_cell),
            'mec_n_hd_95'      : len(hdcells95 & cells_session_mec & parameter_dict_cell),
            'mec_n_border_95'  : len(bordercells95 & cells_session_mec & parameter_dict_cell),
            'mec_n_bv_95'      : len(bvcells95 & cells_session_mec & parameter_dict_cell),
            'mec_n_grid_99'    : len(gridcells99 & cells_session_mec & parameter_dict_cell),
            'mec_n_info_99'    : len(infocells99 & cells_session_mec & parameter_dict_cell),
            'mec_n_hd_99'      : len(hdcells99 & cells_session_mec & parameter_dict_cell),
            'mec_n_border_99'  : len(bordercells99 & cells_session_mec & parameter_dict_cell),
            'mec_n_bv_99'      : len(bvcells99 & cells_session_mec & parameter_dict_cell),
        }

        ##### PAS #######################################################################################################################
        entry_dict_pas = {
            **key,
            'pas_n_cells'      : len(cells_session_pas & parameter_dict_cell),
            'pas_ovc'          : len(ovcells & cells_session_pas & parameter_dict_cell),
            'pas_n_grid_95'    : len(gridcells95 & cells_session_pas & parameter_dict_cell),
            'pas_n_info_95'    : len(infocells95 & cells_session_pas & parameter_dict_cell),
            'pas_n_hd_95'      : len(hdcells95 & cells_session_pas & parameter_dict_cell),
            'pas_n_border_95'  : len(bordercells95 & cells_session_pas & parameter_dict_cell),
            'pas_n_bv_95'      : len(bvcells95 & cells_session_pas & parameter_dict_cell),
            'pas_n_grid_99'    : len(gridcells99 & cells_session_pas & parameter_dict_cell),
            'pas_n_info_99'    : len(infocells99 & cells_session_pas & parameter_dict_cell),
            'pas_n_hd_99'      : len(hdcells99 & cells_session_pas & parameter_dict_cell),
            'pas_n_border_99'  : len(bordercells99 & cells_session_pas & parameter_dict_cell),
            'pas_n_bv_99'      : len(bvcells99 & cells_session_pas & parameter_dict_cell),
        }

        # Set non-OV session entries to None for _ovc 
        # Check whether this session could have OVC entries
        if not (OVC.proj(is_ovc='is_ovc', session_name='base_session') & key):
            entry_dict_mec['mec_ovc'] = None
            entry_dict_pas['pas_ovc'] = None


        # Insert into table 
        self.insert1({**key, **{'n_cells' : len(cell_keys_session)}})
        self.MEC.insert1(entry_dict_mec, ignore_extra_fields=True)
        self.PAS.insert1(entry_dict_pas, ignore_extra_fields=True)

        return 