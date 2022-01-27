import datajoint as dj

from datajoint_utilities.dj_data_copy import db_migration


schema_prefix_update_mapper = {
    'prod_mlims_data': 'group_shared_topopaper_mlims',
    'group_imaging_1b': 'group_shared_topopaper_main_imaging',
    'user_horsto_imaging': 'group_shared_topopaper_horst_imaging',
    'user_horsto_borderscore': 'group_shared_topopaper_borderscore'
    }

all_animals = [
    '90222', '90218', '90647',
    '82913', '88592', '89622',
    '87244', '89841', '60480',
    '87245', '87187', '88106',
    '94557', '97045', '97046',
]

table_block_list = {
    'group_shared_topopaper_mlims': [],
    'group_shared_topopaper_main_imaging': [],
    'group_shared_topopaper_horst_imaging': ['AlignmentPoints',
                                             'AlignmentFOV',
                                             'ScoremapFOV',
                                             'ScoremapFOVMoran',
                                             'ScoremapCorr'],
    'group_shared_topopaper_borderscore': []
}


def main():
    vm_mlims = dj.create_virtual_module('mlims', 'prod_mlims_data')
    vm_imaging = dj.create_virtual_module('imaging', 'group_imaging_1b')

    sessions = (vm_imaging.Session
                @ vm_mlims.Animal & [{'animal_name': n}
                                     for n in all_animals]).fetch('KEY')
    sessions = sessions[::-1]

    batch_size = 5
    for orig_schema_name, cloned_schema_name in schema_prefix_update_mapper.items():
        orig_schema = dj.create_virtual_module(orig_schema_name, orig_schema_name)
        cloned_schema = dj.create_virtual_module(cloned_schema_name, cloned_schema_name)
        # do data copy in batch
        for i in range(0, len(sessions), batch_size):
            db_migration.migrate_schema(orig_schema, cloned_schema,
                                        restriction=sessions[i:i + batch_size],
                                        table_block_list=table_block_list[cloned_schema_name],
                                        allow_missing_destination_tables=True,
                                        force_fetch=True)

    # For tables with renamed parents
    orig_schema_name = 'group_shared_topopaper_horst_imaging'
    cloned_schema_name = schema_prefix_update_mapper[orig_schema_name]
    orig_schema = dj.create_virtual_module(orig_schema_name, orig_schema_name)
    cloned_schema = dj.create_virtual_module(cloned_schema_name, cloned_schema_name)

    for table_name in ['AlignmentPoints', 'AlignmentFOV',
                       'ScoremapFOV', 'ScoremapFOVMoran', 'ScoremapCorr']:
        orig_table = getattr(orig_schema, table_name)
        cloned_table = getattr(cloned_schema, table_name)

        for k in orig_table.fetch('KEY'):
            v = (orig_table & k).fetch1()
            try:
                cloned_table.insert1(
                    v, allow_direct_insert=True, skip_duplicates=True)
            except dj.errors.DataJointError as e:
                print(f'Skipping error: {str(e)}')


if __name__ == '__main__':
    main()
