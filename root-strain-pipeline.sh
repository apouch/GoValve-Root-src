#!/bin/bash

WDIR=$1
FNIMG=$2
FNSEG=$3
NREF=$4
TAG=$5

PATH_GREEDY=/Users/alison/greedy-itk5/BUILD/greedy
PATH_VTKLEVELSET=/Users/alison/cmrep-czi54/build
PATH_MATLAB=/Applications/MATLAB_R2021a.app/bin

# convert segmentation to vtk with label data
FNVTK=$WDIR/segref.vtk
$PATH_VTKLEVELSET/vtklevelset -pl $FNSEG $FNVTK 1

# create medial and boundary meshes of reference segmentation
FNMEDOUT=$WDIR/segref.med.vtk
FNBNDOUT=$WDIR/segref.bnd.vtk
$PATH_MATLAB/matlab -batch "medial_mesh('$FNVTK','$FNMEDOUT','$FNBNDOUT')"

# propagate reference mesh to other frames in the series
python run_propagation_root.py $WDIR $FNIMG $FNSEG $NREF $TAG $PATH_GREEDY $PATH_VTKLEVELSET

# compute root strain
python compute_strain.py ${WDIR}/output/mesh 1 1 2 $NREF

# create 4D root segmentation
python stack_img3d.py ${WDIR}/output $FNIMG $FNSEG $NREF 
