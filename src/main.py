import argparse
from pathlib import Path
from myLog import Log
from preprocessing import PreProcessing
from datetime import timedelta
import time
from clustering_exact import ClusteringExact
from clustering_similar import ClusteringSimilar
from clustering_convergence import ClusteringConvergence
from light_chain import LightChainAnalysis


def main(args):
    """
    CDR3 project: Analyze the frequency of CDR3 combination among various of cancer types from GDC database
    Level 1: Find exact match of CDR3 sequences
    Level 2: Find match of CDR3 sequences having identity of 90%
    Level 3: Refining clusters from level 2 that share the same amino acid group

    INPUT: 
      {cancer_type}
      -i airr.tsv directory path
      -g manifest.txt path (MUST HAVE sample_id and patient_id column)
      -o output directory path
      -l (Optional) log path
      -f (Optional) force to re-run
    """
    start = time.time()

    indir_path = Path(args.indir)
    aa_group_path = Path(args.aminoacid) if args.aminoacid else Path.cwd(
    ).joinpath("reference").joinpath("BLOSUM62.csv")
    if not aa_group_path.is_file():
        raise Exception(
            "NOT FOUND: Amino acid table. Please provide the path to the amino acid table (csv format)")

    gdc_path = Path(args.gdc)

    ##############################################################################################################

    outdir_path = Path(
        args.outdir) if args.outdir else Path.cwd().joinpath("results").joinpath(args.cancer_type)
    outdir_heavy_path = outdir_path.joinpath("heavy")
    outdir_light_path = outdir_path.joinpath("light")

    ##############################################################################################################

    log_path = Path(args.log) if args.log else Path.cwd().joinpath(
        "log").joinpath(args.cancer_type)
    log_path.mkdir(parents=True, exist_ok=True)
    logging = Log(path=log_path)

    logging.info(f"SETTINGS")
    logging.info(f"Input directory path is {indir_path}")
    logging.info(f"Cancer type is {args.cancer_type}")
    logging.info(f"Amino acid score table path is {aa_group_path}")
    logging.info(f"GDC metadata path is {gdc_path}")
    logging.info(f"Log path is {log_path}")
    logging.info(f"Similarity threshold: {args.similarity}")
    logging.info(f"Minimum number of patient in a cluster: {args.patient_min}")
    logging.info(
        f"Output directory is {outdir_path}\n")

    ##############################################################################################################

    try:
        data_heavy_clean_path, data_light_clean_path = PreProcessing(
            indir_path, gdc_path, (outdir_heavy_path, outdir_light_path)).concatenate_airr()

        # heavy chain analysis
        logging.info("HEAVY CHAIN ANALYSIS")
        ClusteringExact(data_heavy_clean_path, args.patient_min,
                        outdir_heavy_path).cluster()
        ClusteringSimilar(
            data_heavy_clean_path, args.patient_min, args.similarity, outdir_heavy_path).cluster()
        group_cluster_path = ClusteringConvergence(data_heavy_clean_path, args.patient_min, aa_group_path, args.similarity,
                                                   outdir_heavy_path).cluster()

        # # light chain analysis
        # logging.info("LIGHT CHAIN ANALYSIS")
        # LightChainAnalysis(group_cluster_path,
        #                    data_heavy_clean_path, data_light_clean_path, args.cancer_type).execute()
    except Exception as e:
        logging.error(e)

    ##############################################################################################################

    running_time = time.time() - start
    logging.info(
        f"Total running time: {str(timedelta(seconds=running_time))}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='CDR3 analysis')
    parser.add_argument('cancer_type', type=str,
                        help='cancer_type')
    parser.add_argument('-i', '--indir',
                        type=str,
                        required=True,
                        help='trust4 output directory (contain _airr.tsv and _report.tsv)')
    parser.add_argument('-g', '--gdc',
                        type=str,
                        required=True,
                        help='path to manifest and metadata gdc file')
    parser.add_argument('-o', '--outdir',
                        type=str,
                        help='output directory')
    parser.add_argument('-l', '--log',
                        type=str,
                        help='output directory')
    parser.add_argument('-a', '--aminoacid',
                        type=str,
                        help='amino acid group mapper')
    parser.add_argument('-p', '--patient_min',
                        type=int,
                        default=2,
                        help='Minimum number of patient in a cluster')
    parser.add_argument('-s', '--similarity',
                        type=float,
                        default=0.8,
                        help='CD-HIT similarity threshold')
    args = parser.parse_args()
    main(args)
