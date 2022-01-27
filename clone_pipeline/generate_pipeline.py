import datajoint as dj
import re
import pathlib

from datajoint_utilities.dj_data_copy import diagram_restriction, db_migration

from .tables_for_bk import tables_for_bk


# ---- reconstruct schema/table code for the selected tables ----
schema_prefix_update_mapper = {
    'prod_mlims_data': 'group_shared_topopaper_mlims',
    'group_imaging_1b': 'group_shared_topopaper_main_imaging',
    'user_horsto_imaging': 'group_shared_topopaper_horst_imaging',
    'user_horsto_borderscore': 'group_shared_topopaper_borderscore'
    }


def main():
    all_tables = []
    for schema_name, selected_table_names in tables_for_bk.items():
        all_tables.extend([f'`{schema_name}`.`{dj.utils.from_camel_case(t).replace(".", "__")}`' for t in selected_table_names])

    schemas_code, schemas_table_definition = diagram_restriction.generate_schemas_definition_code(
        all_tables, schema_prefix_update_mapper,
        verbose=True, save_dir=None)

    # ---- replace blob@store with longblob ----
    for k, code in schemas_code.items():
        schemas_code[k] = re.sub('(blob@\w+)\s', 'longblob', code)

    # ---- write code into .py files ----
    save_dir = pathlib.Path('./clone_pipeline/pipeline_code').absolute()
    save_dir.mkdir(exist_ok=True, parents=True)
    for cloned_schema_name, schema_definition_str in schemas_code.items():
        with open(save_dir / f'{cloned_schema_name}.py', 'wt') as f:
            f.write(schema_definition_str)

    # ---- execute the pipeline codes ----
    # this is where you review the .py files, and run them manually once
    # to get the cloned pipeline initiated
    # (we can code this up, but this should be a manual step with careful inspection)


if __name__ == '__main__':
    main()
