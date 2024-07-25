import os
import pandas as pd

QIITA_SAMPLE_NAME_KEY = 'sample_name'
TAB_SEP = "\t"
SAMPLE_KEY = "Sample"
SAMPLE_NAME_KEY = "Sample Name"
PROJECT_PLATE_KEY = "Project Plate"
SAMPLE_TYPE_KEY = "Sample Type"
PRIMARY_STUDY_KEY = "Primary Qiita Study"
SECONDARY_STUDIES_KEY = "Secondary Qiita Studies"
BLANK_ROOT = "BLANK"
BLANK_SAMPLE_TYPE = "control blank"
STUDY_ID_DELIMITER = ";"
QIITA_STUDY_ID = 'qiita_study_id'
PROJECT_ABBREV_KEY = 'Project Abbreviation'

EXPERIMENT_NAME_KEY = 'Experiment Name'
BIOINFO_KEY = 'Bioinformatics'
CONTACT_KEY = 'Contact'
SAMPLE_CONTEXT_KEY = 'Sample Context'
PROJECT_NAME_KEY = 'Project Name'
HUMAN_FILTERING_KEY = 'HumanFiltering'
EXPT_DESC_KEY = 'experiment_design_description'
CONTAINS_REPS_KEY = 'contains_replicates'
SAMPLE_PROJECT_KEY = 'Sample_Project'
EMAIL_KEY = 'Email'


def get_studies_attr_list(studies_dict, desired_key):
    return [x[desired_key] for x in studies_dict]


def join_dfs_from_files(input_fps, req_cols_to_extract,
                        opt_cols_to_extract=None, unique_cols=None,
                        dtype=None, sep=TAB_SEP):
    if opt_cols_to_extract is None:
        opt_cols_to_extract = []

    if unique_cols is None:
        unique_cols = req_cols_to_extract

    # Per Daniel M. 20240614
    if dtype is None:
        dtype = str

    uniques_required = [x in req_cols_to_extract for x in unique_cols]
    if not all(uniques_required):
        raise ValueError("All unique_cols must be in req_cols_to_extract")

    growing_df = None
    for curr_fp in input_fps:
        if not os.path.isfile(curr_fp):
            raise ValueError(
                "Problem! %s is not a path to a valid file" % curr_fp)

        # load the file into a dataframe, then extract the columns of interest
        curr_df = pd.read_csv(curr_fp, sep=sep, dtype=dtype)
        curr_opt_cols = \
            [x for x in opt_cols_to_extract if x in curr_df.columns]
        curr_cols_to_extract = req_cols_to_extract + curr_opt_cols
        curr_extracted_df = curr_df[curr_cols_to_extract]

        if growing_df is None:
            growing_df = curr_extracted_df
        else:
            not_in_old = \
                set(curr_extracted_df.columns) - set(growing_df.columns)
            not_in_new = \
                set(growing_df.columns) - set(curr_extracted_df.columns)
            for curr_missing in not_in_old:
                growing_df[curr_missing] = None
            for curr_missing in not_in_new:
                curr_extracted_df[curr_missing] = None
            growing_df = pd.concat([growing_df, curr_extracted_df],
                                   ignore_index=True)
        # end if this isn't the first file
    # next file path

    # if any value in the unique_col is duplicated, raise an error
    errs = []
    for curr_unique_col in unique_cols:
        if growing_df[curr_unique_col].duplicated().any():
            dupes = growing_df[curr_unique_col][
                growing_df[curr_unique_col].duplicated()].unique().tolist()
            errs.append(f'Duplicate {curr_unique_col} found in files: {dupes}')
    # next unique_col
    if len(errs) > 0:
        raise ValueError('\n'.join(errs))

    return growing_df


def extend_sample_accession_df(sample_accession_df, studies_info, metadata_df):
    # sample_accession_df should have a 'sample_name' col
    # metadata_df should have 'sample_name' and 'qiita_study_id' cols
    # each entry in studies_info dict should have a key named 'Project Name'
    # and one named 'Project Abbreviation'

    # extract qiita study ids from the Project Name in studies_info entries
    local_studies_info = studies_info.copy()
    for curr_study in local_studies_info:
        curr_project_name = curr_study[PROJECT_NAME_KEY]
        curr_qiita_id = _get_qiita_id_from_sample_project(curr_project_name)
        curr_study[QIITA_STUDY_ID] = curr_qiita_id
    studies_df = pd.DataFrame(local_studies_info)

    # check for qiita_study_ids in studies_df that aren't in the metadata_df
    _check_for_missing_df_ids(studies_df, metadata_df, QIITA_STUDY_ID,
                              'studies', 'metadata')
    # merge the metadata_df with the studies_df on the qiita_study_id
    metadata_plus_df = pd.merge(metadata_df, studies_df, on=QIITA_STUDY_ID)
    # pull the qiita study id off the sample name
    metadata_plus_df[QIITA_SAMPLE_NAME_KEY] = metadata_plus_df.apply(
        lambda x: x[QIITA_SAMPLE_NAME_KEY].replace(f"{x[QIITA_STUDY_ID]}.", ""),
        axis=1
    )

    # I'm pretty sure that actually 'Project Abbreviation' isn't used in
    # *shotgun* processing, but I'm not 100% sure, so I'm keeping it for now.
    # AFAICT, it is used for 16S, to name blanks and katharoseq controls ...
    # in that case it should really be called Plate Abbreviation if more than
    # one project is allowed on a plate ...
    extension_cols = [QIITA_SAMPLE_NAME_KEY, PROJECT_NAME_KEY, PROJECT_ABBREV_KEY]

    # now add the project name and abbreviation to the sample_accession_df
    # by merging on the qiita_id.sample_name
    _check_for_missing_df_ids(
        sample_accession_df, metadata_plus_df, QIITA_SAMPLE_NAME_KEY,
        'sample accession', 'metadata')
    sample_accession_plus_df = pd.merge(
        sample_accession_df, metadata_plus_df[extension_cols],
        on=QIITA_SAMPLE_NAME_KEY)
    return sample_accession_plus_df


def extend_compression_layout_info(compression_layout, studies_info):
    extended_compression_layout = compression_layout.copy()
    # for each dict in compression_layout
    for curr_plate in extended_compression_layout:
        # get its Project Name
        curr_project_name = curr_plate[PROJECT_NAME_KEY]
        # find that Project Name in studies_info
        found_study = None
        for curr_study in studies_info:
            if curr_study[PROJECT_NAME_KEY] == curr_project_name:
                found_study = curr_study
                break
        # next study
        # add the Project Abbreviation to the dict
        if found_study is None:
            raise ValueError(f"{PROJECT_NAME_KEY} '{curr_project_name}' "
                             f"not found in studies_info")
        curr_plate[PROJECT_ABBREV_KEY] = found_study[PROJECT_ABBREV_KEY]
    # next plate

    return extended_compression_layout


def generate_sections_dict(plate_df, studies_info, expt_name, expt_type,
                           expt_version, bioinfo_section_base):
    sections_dict = {
        EXPERIMENT_NAME_KEY: expt_name,
        'SheetType': expt_type,
        'SheetVersion': expt_version,
        'Assay': 'Metagenomic'
    }

    bioinfo_dicts = []
    contacts_dicts = []
    for curr_study_dict in studies_info:
        # bioinformatics dict
        curr_bioinfo_dict = bioinfo_section_base.copy()
        curr_proj_name = curr_study_dict[PROJECT_NAME_KEY]
        curr_qiita_id = _get_qiita_id_from_sample_project(curr_proj_name)
        curr_bioinfo_adds = {
            SAMPLE_PROJECT_KEY: curr_proj_name,
            'QiitaID': curr_qiita_id,
            HUMAN_FILTERING_KEY: curr_study_dict[HUMAN_FILTERING_KEY],
            EXPT_DESC_KEY: curr_study_dict[EXPT_DESC_KEY],
            # Per Charlie 20240715: it is accurate that 'contains_replicates'
            # is added individually to each study dict but has the same value
            # in each because it is really a 384-well-plate-level value.
            # Someday maybe it can move into the header.
            CONTAINS_REPS_KEY: plate_df[CONTAINS_REPS_KEY].all(),
        }
        curr_bioinfo_dict.update(curr_bioinfo_adds)
        bioinfo_dicts.append(curr_bioinfo_dict)

        curr_contact_dict = {
            SAMPLE_PROJECT_KEY: curr_proj_name,
            EMAIL_KEY: curr_study_dict[EMAIL_KEY]
        }
        contacts_dicts.append(curr_contact_dict)
    # next study in the run

    sections_dict[BIOINFO_KEY] = bioinfo_dicts
    sections_dict[CONTACT_KEY] = contacts_dicts
    sections_dict[SAMPLE_CONTEXT_KEY] = _generate_sample_context(plate_df)

    return sections_dict


# Generate the metadata dictionary SampleContext section, which looks for
# example like the below, and add it to the metadata dictionary
#
#    'SampleContext': [
#        {
#            'Sample_Name': 'BLANK.NPH.4.G11',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '12000;10981',
#            'Sample_Type': 'control blank'
#        },
#        {
#            'Sample_Name': 'BLANK.NPH.8.G10',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '',
#            'Sample_Type': 'control blank'
#        },
#        {
#            'Sample_Name': 'BLANK.NPH.18.G01',
#            'PrimaryQiitaStudy': '12986',
#            'SecondaryQiitaStudies': '10981',
#            'Sample_Type': 'control blank'
#        },
#    ]
def _generate_sample_context(the_plate_df, blanks_mask=None):
    if blanks_mask is None:
        blanks_mask = _get_blanks_mask_by_name_format(the_plate_df)

    # get all the unique Project Plate values for the blanks
    names_of_plates_w_blanks = \
        the_plate_df.loc[blanks_mask, PROJECT_PLATE_KEY].unique()

    # start building the SampleContext object
    blanks_context_df = the_plate_df.loc[
        blanks_mask,
        [SAMPLE_KEY, PROJECT_NAME_KEY, PROJECT_PLATE_KEY]].copy()
    blanks_context_df.rename(columns={SAMPLE_KEY: SAMPLE_NAME_KEY},
                             inplace=True)
    blanks_context_df[PRIMARY_STUDY_KEY] = \
        blanks_context_df.apply(_get_primary_study, axis=1)
    blanks_context_df[SAMPLE_TYPE_KEY] = BLANK_SAMPLE_TYPE

    # for each plate in names_of_plates_w_blanks
    for name_of_curr_plate_w_blanks in names_of_plates_w_blanks:
        # note that this mask is for ALL of plate_df, not just blank rows
        curr_plate_mask = \
            the_plate_df[PROJECT_PLATE_KEY] == name_of_curr_plate_w_blanks
        # get all the unique Project Name values for plate_df rows that have
        # Project Plate value == name_of_curr_plate_w_blanks
        projects_on_curr_plate = the_plate_df.loc[
            curr_plate_mask, PROJECT_NAME_KEY].unique().tolist()
        # add projects_on_curr_plate to blanks_context_df for
        # this plate's blanks
        curr_plate_blanks_mask = \
            blanks_context_df[PROJECT_PLATE_KEY] == name_of_curr_plate_w_blanks

        # for each row in blanks_context_df that has the current plate name,
        # remove that row's project name from projects_on_curr_plate list to
        # generate the list of secondary studies for that blank
        blanks_context_df.loc[
            curr_plate_blanks_mask, SECONDARY_STUDIES_KEY] = (
            blanks_context_df.apply(
                _get_secondary_studies,
                all_projects_on_plate=projects_on_curr_plate,
                axis=1)
        )
    # next plate name of plate w blanks

    # clean up by removing the Project Plate and Project Name columns
    blanks_context_df.drop(columns=[PROJECT_PLATE_KEY, PROJECT_NAME_KEY],
                           inplace=True)

    # turn blanks_context_df into a list of dictionaries
    blanks_context_list = blanks_context_df.to_dict(orient="records")
    return blanks_context_list


def _get_blanks_mask_by_name_format(a_plate_df):
    return a_plate_df[SAMPLE_KEY].str.startswith(BLANK_ROOT) == True


def _get_qiita_id_from_sample_project(sample_project_str):
    return sample_project_str.split("_")[-1]


def _check_for_missing_df_ids(subset_df, superset_df, id_col,
                              sub_name, sup_name):
    # if there are any qiita_study_ids in the studies_df that aren't in the
    # metadata_df, raise an error
    found_mask = subset_df[id_col].isin(superset_df[id_col])
    if not found_mask.all():
        missing_vals = subset_df.loc[~found_mask, id_col].unique()
        raise ValueError(f'Some {id_col} values in the {sub_name} dataframe '
                         f'are not in the {sup_name} dataframe: '
                         f'{", ".join(missing_vals)}')


def _get_primary_study(curr_row):
    return _get_qiita_id_from_sample_project(curr_row[PROJECT_NAME_KEY])


def _get_secondary_studies(curr_row, all_projects_on_plate):
    row_project = curr_row[PROJECT_NAME_KEY]

    # NB: below will error if row_project is not in
    # projects_on_curr_plate; this is as it should be because
    # that situation should never arise and if it does, we need to
    # stop and figure out why
    secondary_projects = all_projects_on_plate.copy()
    secondary_projects.remove(row_project)
    secondary_projects = sorted(secondary_projects)

    # now split each secondary_projects on _ and make a list of the
    # last element of each split
    secondary_qiita_ids = [_get_qiita_id_from_sample_project(x)
                           for x in secondary_projects]

    secondary_qiita_ids_str = STUDY_ID_DELIMITER.join(secondary_qiita_ids)
    return secondary_qiita_ids_str
