python run_regression.py --feature_list presurgical
python run_regression.py --feature_list presurgical_reduced
python run_regression.py --feature_list postsurgical

sbatch hpc/run_array_cpu.sbatch --array 0-999 --export=feature_list="presurgical"
sbatch hpc/run_array_cpu.sbatch --array 0-999 --export=feature_list="presurgical_reduced"
sbatch hpc/run_array_cpu.sbatch --array 0-999 --export=feature_list="postsurgical"

