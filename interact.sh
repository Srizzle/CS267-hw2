#!/bin/bash
#SBATCH -A cc3uv3p # 2017 XSEDE Applications of Parallel Computing Course Allocation
#SBATCH -J auto-particle-gpu
#SBATCH -o auto-particle-gpu.stdout
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:k80:1
#SBATCH -t 00:10:00
#SBATCH -N 1         
interact -A cc3uv3p -p GPU-shared --gres=gpu:k80:1 -t 01:00:00 -N 1