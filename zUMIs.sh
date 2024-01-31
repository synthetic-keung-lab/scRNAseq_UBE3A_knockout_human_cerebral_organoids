#!/bin/tcsh
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -x
#BSUB -R "rusage[mem=80GB]"
#BSUB -W 5760
limit descriptors 8191
/usr/local/usrapps/zUMIs/v2.9.4d/zUMIs/zUMIs.sh -c -y zUMIs.yaml #execute zUMIs
