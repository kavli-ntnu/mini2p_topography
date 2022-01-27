import datajoint as dj


# -------------- group_shared_topopaper_mlims --------------


schema = dj.Schema('group_shared_topopaper_mlims')



@schema
class AnimalColor(dj.Lookup):
    definition = """
    color                : varchar(16)                  
    """


@schema
class Strain(dj.Lookup):
    definition = """
    strain               : varchar(64)                  
    """


@schema
class Gene(dj.Lookup):
    definition = """
    gene                 : varchar(64)                  # gene name
    ---
    gene_description     : varchar(4096)                
    """


@schema
class Animal(dj.Manual):
    definition = """
    animal_id            : char(16)                     
    datasource_id=0      : int                          
    ---
    animal_species       : varchar(16)                  
    animal_name          : varchar(128)                 
    animal_sex           : enum('M','F','U')            
    animal_dob=null      : date                         
    -> AnimalColor
    animal_notes=null    : varchar(4096)                
    INDEX (animal_name, animal_species)
    """

    class BackgroundStrain(dj.Part):
        definition = """
        -> Animal
        ---
        background_strain    : varchar(128)                 # background strain (e.g. cross)
        """

    class Generation(dj.Part):
        definition = """
        # animal generation information
        -> Animal
        ---
        generation           : varchar(64)                  # generation information
        """

    class Litter(dj.Part):
        definition = """
        # litter information
        -> Animal
        ---
        litter               : varchar(64)                  # animal litter
        """

    class Purpose(dj.Part):
        definition = """
        # purpose of the animal
        -> Animal
        ---
        purpose              : varchar(64)                  # animal purpose
        """

    class Strain(dj.Part):
        definition = """
        -> Animal
        ---
        -> Strain
        """


@schema
class Genotype(dj.Manual):
    definition = """
    -> Animal
    -> Gene
    ---
    genotype             : varchar(128)                 # genotype constitution
    """


