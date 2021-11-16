#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 11:49:45 2021

@author: alison
"""

import sys
import SimpleITK as sitk

fnimg = sys.argv[1]

# 4D reference image
img_ref = sitk.ReadImage(fnimg)

# Print frame time
ft = img_ref.GetSpacing()[3]
print(ft)