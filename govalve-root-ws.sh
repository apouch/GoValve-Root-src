#!/bin/bash

WS=$1
WDIR=$2

PATH_SNAP="/home/apouch/src/itksnap-4.0.0-alpha-20211103-Linux-gcc64/bin"

fn_img4D=$(${PATH_SNAP}/itksnap-wt -P -i $WS -llf "Source")
fn_seg4D=$(${PATH_SNAP}/itksnap-wt -P -i $WS -llf "Segmentation")
f_ref=$(${PATH_SNAP}/itksnap-wt -P -i $WS -timepoints-pick-by-tag "Reference Frame")
f_start=$(${PATH_SNAP}/itksnap-wt -P -i $WS -timepoints-pick-by-tag "Start Cardiac Cycle")
f_open=$(${PATH_SNAP}/itksnap-wt -P -i $WS -timepoints-pick-by-tag "Valve Opening")
f_close=$(${PATH_SNAP}/itksnap-wt -P -i $WS -timepoints-pick-by-tag "Valve Closing")
f_end=$(${PATH_SNAP}/itksnap-wt -P -i $WS -timepoints-pick-by-tag "End Cardiac Cycle")
tag="root"

# Extract 3D image
fn_seg3D=$WDIR/seg_ref3D.nii.gz
c4d $fn_seg4D -slice w $(($f_ref-1)) -o $fn_seg3D

bash /home/apouch/GoValve-Root-src/root-strain-pipeline.sh $WDIR $fn_img4D $fn_seg3D $f_ref $f_start $f_open $f_close $f_end $tag
zip -r $WDIR/output/mesh/output-meshes.zip $WDIR/output/mesh/Strains 

RESULT_WSP=$WDIR/result.itksnap
#${PATH_SNAP}/itksnap-wt -layers-set-main $fn_img4D \
#           -layers-set-seg $WDIR/output/seg4d.nii.gz \
#           -o $RESULT_WSP

${PATH_SNAP}/itksnap-wt -i $WS \
	-layers-set-main $fn_img4D \
	-layers-set-seg $WDIR/output/seg4d.nii.gz \
	-o $RESULT_WSP
