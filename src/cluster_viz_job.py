import subprocess
from pathlib import Path

log_dir = Path("/rsrch4/home/mol_cgenesis/nkdang/CDR3/log/viz")
log_dir.mkdir(parents=True, exist_ok=True)

path_lsf = Path("/rsrch4/home/mol_cgenesis/nkdang/CDR3/lsf/viz")
path_lsf.mkdir(parents=True, exist_ok=True)

for dir in Path("/rsrch4/home/mol_cgenesis/EMC_BIC_rsrch4/nkdang/CDR3/results").glob("*/heavy/*"):
    if dir.is_dir():
        name = "_".join(str(dir).split("/")[-3:])
        command = "#BSUB -W 12:00\n" \
            "#BSUB -q medium\n" \
            f"#BSUB -n 4\n" \
            f"#BSUB -M 50\n" \
            f"#BSUB -R rusage[mem=50]\n" \
            f"#BSUB -o {log_dir.joinpath(name)}\n" \
            "#BSUB –u nkdang@mdanderson.org\n" \
            f"#BSUB -J VIZ_{dir}\n" \
            f"python /rsrch4/home/mol_cgenesis/nkdang/CDR3/src/cluster_viz.py {dir}\n" \

        with open(path_lsf.joinpath(f"{name}.lsf"), "w+") as out:
            out.writelines(command)
        subprocess.Popen(
            f"bsub < {path_lsf.joinpath(f'{name}.lsf')}", shell=True)
