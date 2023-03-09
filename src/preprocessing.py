import pandas as pd
from pathlib import Path
import logging
from p_tqdm import p_map

counter = 0  # initial value


def func(row):
    global counter
    counter += 1
    return counter


class PreProcessing():
    """
    Merge all the _airr.tsv files in the input directory into a single _airr.tsv file
    Read manifest file (file must contain patient_id column and sample_id)
    Add some additional columns: v_call_ori, j_call_ori, junction_aa_length, patient_id, sample_id, cancer_type
    Drop columns: sequence_alignment, germline_alignment, cell_id
    """

    def __init__(self, indir_path: Path, gdc_path: Path, outdir_path: tuple[Path, Path]) -> None:
        logging.info("PREPROCESSING ...")
        self.indir_path = indir_path

        # create mapper from case_id to patient_id
        manifest = pd.read_csv(gdc_path)

        manifest = manifest.sort_values(
            ['patient_id'])

        self.mapper_id_patient = dict(
            zip(manifest['sample_id'].astype(str), manifest['patient_id'].astype(str)))

        self.outdir_heavy_path, self.outdir_light_path = outdir_path
        self.outdir_heavy_path.mkdir(parents=True, exist_ok=True)
        self.outdir_light_path.mkdir(parents=True, exist_ok=True)

    def __read_airr__(self, filename: Path) -> pd.DataFrame:
        """
        converts a filename to a pandas dataframe
        """
        if not filename.parent.name in self.mapper_id_patient.keys():
            return
        patient_id = self.mapper_id_patient[filename.parent.name]

        dtypes = {'sequence_id': 'str',
                  'sequence': 'str',
                  'rev_comp': 'str',
                  'productive': 'str',
                  'v_call': 'str',
                  'd_call': 'str',
                  'j_call': 'str',
                  'c_call': 'str',
                  'cdr1': 'str',
                  'cdr2': 'str',
                  'junction': 'str',
                  'junction_aa': 'str',
                  'v_cigar': 'str',
                  'd_cigar': 'str',
                  'j_cigar': 'str',
                  'v_identity': 'str',
                  'j_identity': 'str',
                  'complete_vdj': 'str',
                  'consensus_count': 'int64',
                  'sample_id': 'str',
                  'cancer_type': 'str'}

        airr = pd.read_csv(filename, sep='\t', dtype=dtypes)
        airr = airr.dropna(
            subset=['v_call', 'j_call', 'junction_aa'], how='any')

        assert airr[['v_call', 'j_call', 'junction_aa']].isna(
        ).sum().sum() == 0, f"{filename} There're NaN values"

        airr['v_call_ori'] = airr['v_call']
        airr['v_call'] = airr['v_call'].apply(lambda x: str(x).split("*")[0])

        airr['j_call_ori'] = airr['j_call']
        airr['j_call'] = airr['j_call'].apply(lambda x: str(x).split("*")[0])

        assert airr[['v_call', 'j_call', 'junction_aa']].isna(
        ).sum().sum() == 0, f"{filename} There're NaN values"

        airr['sample_index'] = airr.index
        airr['sample_index'] = airr['sample_index'].astype(str)
        airr['sample_id'] = filename.parent.name
        # airr['normalized_count'] = airr['consensus_count'] / \
        #     airr['consensus_count'].sum()*1000

        airr['IG_subtypes'] = airr['c_call'].str[:4]
        airr = airr[['sample_id', 'sample_index', 'v_call',
                     'j_call', 'junction_aa', 'IG_subtypes', 'consensus_count']]
        airr['patient_id'] = patient_id

        return airr

    def concatenate_airr(self) -> None:
        """
        combine all _airr.tsv files into a single file
        """
        # return self.outdir_heavy_path.joinpath("data_clean.feather"), self.outdir_light_path.joinpath("data_clean.feather")
        outdir_others_path = self.outdir_heavy_path.parent.joinpath("others")
        outdir_others_path.mkdir(parents=True, exist_ok=True)
        if not (self.outdir_heavy_path.joinpath("data_clean.csv").is_file() and self.outdir_light_path.joinpath("data_clean.csv").is_file() and outdir_others_path.joinpath("data_clean.csv").is_file()) or True:
            logging.info(
                f"Merge all _airr.tsv files from {self.indir_path}")
            total_cases = len(list(self.indir_path.glob('**/*_airr.tsv')))
            logging.info(f"Total sample files: {total_cases:,d}")
            logging.info(
                f"Total unique cancer patients: {len(self.mapper_id_patient):,d}")

            # have your pool map the file names to dataframes
            df_list = p_map(self.__read_airr__, self.indir_path.glob(
                "**/*_airr.tsv"), total=total_cases)

            # reduce the list of dataframes to a single dataframe
            df_combined = pd.concat(df_list, ignore_index=True)

            logging.info(f"Total records: {len(df_combined):,d}")

            df_heavy = df_combined[df_combined['v_call'].str.startswith('IGH')]
            df_light = df_combined[(df_combined['v_call'].str.startswith(
                'IGK')) | (df_combined['v_call'].str.startswith('IGL'))]
            df_others = df_combined[~(df_combined.index.isin(
                df_heavy.index) | df_combined.index.isin(df_light.index))]

            logging.info(f"Writing heavy chain data IGH ...")
            df_heavy = df_heavy.reset_index(drop=True)
            df_heavy_no_dup = df_heavy.groupby(['sample_id', 'v_call', 'j_call', 'junction_aa', 'patient_id', 'IG_subtypes']).aggregate(
                {"consensus_count": "sum", 'sample_index': "_".join}).reset_index()
            df_heavy_no_dup.to_csv(self.outdir_heavy_path.joinpath(
                "data_clean.csv"), index=False)
            df_heavy_no_dup.to_feather(self.outdir_heavy_path.joinpath(
                "data_clean.feather"))
            logging.info(
                f"Heavy chain records: {len(df_heavy):,d} | {len(df_heavy)/len(df_combined)*100:.2f}%")
            logging.info(
                f"After aggregate duplicated records: {len(df_heavy_no_dup):,d}")

            logging.info(f"Writing light chain data IGK and IGL ...")
            df_light = df_light.reset_index(drop=True)
            df_light_no_dup = df_light.groupby(['sample_id', 'v_call', 'j_call', 'junction_aa', 'patient_id', 'IG_subtypes']).aggregate(
                {"consensus_count": "sum", 'sample_index': "_".join}).reset_index()
            df_light_no_dup.to_csv(self.outdir_light_path.joinpath(
                "data_clean.csv"), index=False)
            df_light_no_dup.to_feather(self.outdir_light_path.joinpath(
                "data_clean.feather"))
            logging.info(
                f"Light chain records: {len(df_light):,d} | {len(df_light)/len(df_combined)*100:.2f}%")
            logging.info(
                f"After aggregate duplicated records: {len(df_light_no_dup):,d}")

            logging.info(f"Writing others data ...")
            df_others = df_others.reset_index(drop=True)
            df_others.to_csv(outdir_others_path.joinpath(
                "data_clean.csv"), index=False)
            df_others.to_feather(outdir_others_path.joinpath(
                "data_clean.feather"))
            logging.info(
                f"Other records: {len(df_others):,d} | {len(df_others)/len(df_combined)*100:.2f}%")

        logging.info(
            f"Finished preprocessing data.\n")

        return self.outdir_heavy_path.joinpath("data_clean.feather"), self.outdir_light_path.joinpath("data_clean.feather")
