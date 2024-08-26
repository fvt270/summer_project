#!/bin/sh
#SBATCH --job-name=hiv_LPV
#SBATCH --time=48:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab_ib

conda activate /groups/sbinlab/software/PRISM_tools/py3_ros_ddG_env

dir_py=/groups/sbinlab/fvt270/PRISM/software/rosetta_ddG_pipeline
dir_run=/groups/sbinlab/fvt270/summer_project/drugs/LPV

### LPV
python3 $dir_py/run_pipeline.py -s $dir_run/1mui_aligned.pdb -o $dir_run/output_LPV -i create -mm mut_file -m $dir_run/LPV_mutfile.txt --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

### Stability_calculations -i proceed instead of -i create
python3 $dir_py/run_pipeline.py -s $dir_run/1mui_aligned.pdb -o $dir_run/output_LPV -i proceed -mm mut_file -m $dir_run/LPV_mutfile.txt --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True

### ONLY ddg_calculation
#python3 $dir_py/run_pipeline.py -s $dir_run/1mui_aligned.pdb -o $dir_run/output_LPV -i ddg_calculation -mm mut_file -m $dir_run/LPV_mutfile.txt --chainid AB --run_struc AB --overwrite_path True --slurm_partition sbinlab_ib --ligand True --dump_pdb True
