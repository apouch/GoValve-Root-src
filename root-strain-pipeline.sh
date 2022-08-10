#!/bin/bash

WDIR=$1
FNIMG=$2
FNSEG=$3
NREF=$4
NSTART=$5
NOPEN=$6
NCLOSE=$7
NDONE=$8
TAG=$9

echo "Working directory: $WDIR"
echo "Image filename: $FNIMG"
echo "Segmentation filename: $FNSEG"
echo "Reference frame number: $NREF"
echo "Start frame number: $NSTART"
echo "Open frame number: $NOPEN"
echo "Close frame number: $NCLOSE"
echo "End frame number: $NDONE"
echo "Tag: $TAG"

PATH_GREEDY=/home/apouch/src/itksnap-4.0.0-alpha-20211103-Linux-gcc64/bin
PATH_VTKLEVELSET=/home/apouch/build/cmrep-dev-czi-54
PATH_MATLAB=/usr/local/bin
PATH_SCRIPTS=/home/apouch/GoValve-Root-src

cd $PATH_SCRIPTS

# get frame time of 4D image
FT=$(python get_frametime.py ${FNIMG})
echo "Frame time is ${FT} ms"

# extract 3D segmentation

# convert segmentation to vtk with label data
FNVTK=$WDIR/segref.vtk
#$PATH_VTKLEVELSET/vtklevelset -pl $FNSEG $FNVTK 1

# create medial and boundary meshes of reference segmentation
FNMEDOUT=$WDIR/segref.med.vtk
FNBNDOUT=$WDIR/segref.bnd.vtk
#$PATH_MATLAB/matlab -batch "medial_mesh('$FNVTK','$FNMEDOUT','$FNBNDOUT')"

# propagate reference mesh to other frames in the series
#python run_propagation_root.py $WDIR $FNIMG $FNSEG $NREF $NSTART $NDONE $TAG $PATH_GREEDY $PATH_VTKLEVELSET

# construct medial surfaces from boundary meshes
MESHDIR=$WDIR/output/mesh
#$PATH_MATLAB/matlab -batch "medial_recon_from_bnd('$MESHDIR','$TAG','$FNBNDOUT','$FNMEDOUT','$NREF','$NSTART','$NDONE')"
#cp $FNMEDOUT $MESHDIR/seg_${TAG}_med_recon_${NREF}.vtk

# compute root strain
python compute_strain.py ${WDIR}/output/mesh $FT $NOPEN $NCLOSE $NREF

# create 4D root segmentation
echo "${WDIR}/output"
echo "Image filename is: $FNIMG"
echo "Segmentation filename is: $FNSEG"
echo "Reference frame is: $NREF"
python stack_img3D.py ${WDIR}/output $FNIMG $FNSEG $NREF 
