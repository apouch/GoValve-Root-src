#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 15:03:33 2021

@author: alison
"""

import SimpleITK as sitk
import os
import glob
import sys

WDIR = sys.argv[1]
FNIMG = sys.argv[2]
FNSEG = sys.argv[3]
NREF = int(sys.argv[4])

print("Stacking segmentation volumes")

# list of propagated segmentations
fnseg_prop = sorted(glob.glob(os.path.join(WDIR,"seg*_to_*reslice.nii.gz")))

# number of frames
img_ref = sitk.ReadImage(FNIMG)
nf = img_ref.GetSize()[3]
print(nf)

# reference segmentation and blank copy
seg_ref = sitk.ReadImage(FNSEG)
seg_ref_blank = 0*seg_ref

# append volumes
vol = []
for i in range(nf):
    j = i + 1
    if j == NREF:
        vol.append(seg_ref)
    else:
        s = "seg*_to_" + str(j) + "_*reslice.nii.gz"
        fn = sorted(glob.glob(os.path.join(WDIR,s)))
        if len(fn) == 0:
            vol.append(seg_ref_blank)
        else:
            vol.append(sitk.ReadImage(fn[0]))
    
# create 4D segmentation series    
seg4d = sitk.JoinSeries(vol)
print(seg4d.GetSize())

# write 4D segmentation image
writer = sitk.ImageFileWriter()
writer.SetFileName(os.path.join(WDIR,"seg4d.nii.gz"))
writer.Execute(seg4d)
