import subprocess
from pathlib import Path

log_dir = Path("/rsrch4/home/mol_cgenesis/nkdang/CDR3/hpc/results")
log_dir.mkdir(parents=True, exist_ok=True)

path_lsf = Path("/rsrch4/scratch/mol_cgenesis/nkdang/CDR3/lsf/results")
path_lsf.mkdir(parents=True, exist_ok=True)

for file in Path("/rsrch4/home/mol_cgenesis/nkdang/CDR3/gdc/portal").glob("*"):
    if file.is_dir():
        cancer = file.name

        out_dir = Path(
            f"/rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results").joinpath(cancer)
        # out_dir = Path(
        #     f"/rsrch4/scratch/mol_cgenesis/nkdang/CDR3/results").joinpath(cancer)
        out_dir.mkdir(parents=True, exist_ok=True)

        command = "#BSUB -W 12:00\n" \
            "#BSUB -q medium\n" \
            f"#BSUB -n 8\n" \
            f"#BSUB -M 100\n" \
            f"#BSUB -R rusage[mem=100]\n" \
            "#BSUB –u nkdang@mdanderson.org\n" \
            f"#BSUB -J main_{cancer}\n" \
            f"#BSUB -o {log_dir.joinpath(cancer)}\n" \
            f"#BSUB –cwd /rsrch4/home/mol_cgenesis/nkdang/CDR3\n" \
            f"python /rsrch4/home/mol_cgenesis/nkdang/CDR3/src/main.py {cancer} -i /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/trust4/portal/{cancer} -g /rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/gdc/portal/{cancer}/gdc_manifest.2022-10-08.txt -o {out_dir}\n" \

        with open(path_lsf.joinpath(f"{cancer}.lsf"), "w+") as out:
            out.writelines(command)
        subprocess.Popen(
            f"bsub < {path_lsf.joinpath(f'{cancer}.lsf')}", shell=True)
