import datajoint as dj

#### LOAD DATABASE #########################################
from .dj_conn import *
imhotte = dj.schema(horst_imaging_db)
 

@imhotte 
class AnatomicalMaskParams(dj.Lookup):
    definition = """
    # LUT for anatomical masks drawn in Napari to identify subregions in FOV
    timestamp_mask_lookup: timestamp  # Auto created timestamp
    ---    
    mec_label            : tinyint   # Medial entorhinal cortex (MEC) label name (integer)
    pas_label = NULL     : tinyint   # Parasubiculum label name (integer)
    rsa_label = NULL     : tinyint   # Retrosplenial agranular cortex (integer)
    prh_label = NULL     : tinyint   # Perirhinal cortex
    """ 


@imhotte 
class AnatomicalMask(dj.Manual):
    definition = """
    # Anatomical mask identifying anatomical regions in FOV
    -> ProjectionCorr
    ---    
    -> AnatomicalMaskParams
    anat_mask     : blob@imgstore     # Anatomical mask for regions in FOV
    """
    
@imhotte
class RoisCorrBrainLoc(dj.Computed):
    definition = """
    # Cell IDs and anatomical location 
    -> AnatomicalMask
    -> RoisCorr
    ---
    """
    class MEC(dj.Part):
        definition = """
        # Cells in MEC
        -> master
        ---
        """
    class PAS(dj.Part):
        definition = """
        # Cells in Parasubiculum
        -> master
        ---
        """
    class RSA(dj.Part):
        definition = """
        # Cells in Retrosplenial / Agranular cortex
        -> master
        ---
        """
    class PRH(dj.Part):
        definition = """
        # Cells in perirhinal cortex
        -> master
        ---
        """
    class Undefined(dj.Part):
        # Has label "0" in anatomical mask 
        definition = """
        # Cells elsewhere (not defined by any anatomical label)
        -> master
        ---
        """        
        
    def make(self,key):
        anatomical_mask = (AnatomicalMask & key).fetch1('anat_mask')
        anatomical_lut  = (AnatomicalMaskParams & key).fetch1()
        
        # Sanity check - all labels present in LUT?
        unique_labels_map = set(anatomical_mask.astype(int).flatten())
        unique_labels_map.discard(0) # Don't count 0 because it is "undefined" (see below)
        unique_labels_lut = set([value for value in anatomical_lut.values() if isinstance(value,int)])
        if not all(elem in unique_labels_lut for elem in unique_labels_map): 
            raise KeyError('Labels found in anatomical map not found in LUT')
        
        self.insert1(key)
        
        cell = (RoisCorr.proj('center_x_corr','center_y_corr') & key).fetch1()
        
        if (cell['center_x_corr'] > anatomical_mask.shape[1]) or (cell['center_y_corr'] > anatomical_mask.shape[0]):
            self.Undefined.insert1({**key,**cell},  ignore_extra_fields=True)

        elif anatomical_mask[int(cell['center_y_corr']),int(cell['center_x_corr'])] == anatomical_lut['mec_label']:
            self.MEC.insert1({**key,**cell},  ignore_extra_fields=True)

        elif anatomical_mask[int(cell['center_y_corr']),int(cell['center_x_corr'])] == anatomical_lut['pas_label']:
            self.PAS.insert1({**key,**cell},  ignore_extra_fields=True)

        elif anatomical_mask[int(cell['center_y_corr']),int(cell['center_x_corr'])] == anatomical_lut['rsa_label']:
            self.RSA.insert1({**key,**cell},  ignore_extra_fields=True)

        elif anatomical_mask[int(cell['center_y_corr']),int(cell['center_x_corr'])] == anatomical_lut['prh_label']:
            self.PRH.insert1({**key,**cell},  ignore_extra_fields=True)

        elif anatomical_mask[int(cell['center_y_corr']),int(cell['center_x_corr'])] == 0: 
            self.Undefined.insert1({**key,**cell},  ignore_extra_fields=True)
        else: 
            raise KeyError(f'Label {anatomical_mask[int(cell["center_y_corr"]),int(cell["center_x_corr"])]} not known')



