#!/bin/bash
#SBATCH --job-name=2015_01_4hr_150117_SMOTEs
#SBATCH --output=/fred/oz100/sgoode/SMOTEs/Scripts/slurm_logs/2015_01_4hr_150117_SMOTEs.out
#SBATCH --error=/fred/oz100/sgoode/SMOTEs/Scripts/slurm_logs/2015_01_4hr_150117_SMOTEs.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=3gb
#SBATCH --partition=skylake
#SBATCH --gres=gpu:1

python SMOTE_pipe_wrapper.py -v 2015_01_4hr_150117.list /fred/oz100/NOAO_archive/archive_NOAO_data/scripts/create_lc/image_mjd_lists/FINAL_2015_01_4hr_150117.ascii /fred/oz100/NOAO_archive/archive_NOAO_data/data_outputs/2015/01/4hr/g_band/single/ /fred/oz100/sgoode/SMOTEs/Data/2015_01_4hr/


