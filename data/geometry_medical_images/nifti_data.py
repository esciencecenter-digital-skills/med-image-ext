#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nifti_data

Outputs data related to NIfTI file


Copyright 2021-2022, David Atkinson, University College London

"""


import nibabel   # not in Anaconda - did a pip install nibabel in console.
import os

# Data location and file names
fpath = os.path.join(os.getcwd() , 'NIFTI')

fn_cor = "OBJECT_phantom_T2W_TSE_Cor_14_1.nii" # Coronal

ffn = os.path.join(fpath, fn_cor)
ns = nibabel.load(ffn)
print(ns.affine)
print(ns.header)

ad = ns.get_fdata()
print(ad.flags)
print(ns.shape)
print(ad.shape)













