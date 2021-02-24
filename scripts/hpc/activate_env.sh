#!/usr/bin/bash
module load miniconda/3
source activate ~/.conda/envs/meld
export FREESURFER_HOME=/usr/local/software/freesurfer/6.0.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh

vglrun freesurfer

