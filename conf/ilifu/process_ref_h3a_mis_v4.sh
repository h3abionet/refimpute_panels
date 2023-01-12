#!/usr/bin/env bash
#SBATCH --partition=Main
#SBATCH --nodes=1 --ntasks=2 --mem=7000
#SBATCH --time=240:00:00
#SBATCH --job-name="ref_h3a_mis"
#SBATCH --mail-user=mbymam001@myuct.ac.za
#SBATCH --mail-type=BEGIN,END,FAIL
#SBTACH -o /cbio/dbs/refpanels/H3AR4x/nextflow.out

cd /cbio/users/mamana/refimpute

nextflow \
    ~/refimpute/process_ref_h3a_mis.nf \
    -c /cbio/projects/001/clients/refimpute/process_ref_h3a_mis_v4.config \
    -profile singularity,slurm \
    -resume


# rsync -arvP -e "ssh -i /users/mamana/.ssh/id_rsa" 1.0.0 ubuntu@154.114.37.238:~/