#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geomdemo
Demonstration to accompany 'Geometry in Medical Imaging: DICOM and NIfTI Formats'
Creates a plot of slice positions and two 3D plots demonstrating geometry 
and coordinates.

Data available from Zenodo:
    https://doi.org/10.5281/zenodo.6466491
    
See similar MATLAB code geomdemo.m

User needs to edit the file path variable or run from downloaded folder. 


Copyright 2021-2022, David Atkinson, University College London

"""

# To interact with undocked 3D plots in Spyder, in the iPython console type 
# %matplotlib qt

import pydicom
import matplotlib.pyplot as plt
import numpy as np
import os

# Data location and file names
fpath = os.path.join(os.getcwd() , 'DICOM')

fn_cor = "IM_0002" # Coronal
fn_tra = "IM_0005" # Transverse
fn_sag = "IM_0008" # Sagittal
fn_obl = "IM_0011" # Oblique

def dslice_info(ax, fn, fpath):
    ffn = os.path.join(fpath, fn)
    
    ds = pydicom.dcmread(ffn)

    dpfg = ds.PerFrameFunctionalGroupsSequence 
    
    nframe = len(dpfg)
    proj = np.zeros((nframe))
        
    # Get IOP from first frame
    IOPMV = dpfg[0].PlaneOrientationSequence[0].ImageOrientationPatient
    IOP = np.array([float(IOPMV[0]), float(IOPMV[1]), float(IOPMV[2]), float(IOPMV[3]), float(IOPMV[4]), float(IOPMV[5]),])
    C = np.cross(IOP[0:3], IOP[3:6])
        
    for frame in range(0,nframe):
        IPPMV = dpfg[frame].PlanePositionSequence[0].ImagePositionPatient
        IPP = np.array([float(IPPMV[0]), float(IPPMV[1]), float(IPPMV[2])])
        # Could read IOP for this frame and check same as first frame (parallel stack)
            
        proj[frame] = np.dot(C, IPP)
    # end for
        
    print('Distance between slice 1 and 2 centres is: ' +  str(proj[1] - proj[0]) + ' mm ')

    if (proj[1] < proj[0]):
        print('WARNING: First and 2nd slices are not in ascending order for: ' + fn)
    
          
    # See MATLAB code for additional suggestions for checking slices parallel and equally spaced
    ax.plot(range(0,nframe), proj, label=fn)
    ax.set_xlabel('Frame number (0-based)')
    ax.set_ylabel('Projection')
    ax.legend()
        
def dslice(ax, fn, fpath, frame, splot):
    # Puts one slice on a 2x2 plot
    ffn = os.path.join(fpath, fn)
    ds = pydicom.dcmread(ffn)
    data = ds.pixel_array
    img = data[frame]
    plt.subplot(2,2,splot)
    plt.imshow(img)
    plt.title(fn)
    

def ddisp(ax, fn, fpath, method, frame):
    # ddisp Dicom display in 3D
    ffn = os.path.join(fpath, fn)
    
    ds = pydicom.dcmread(ffn)

    dpfg = ds.PerFrameFunctionalGroupsSequence 
  

    # pydicom returns MultiValue types that need converting below to float arrays
    PSMV  = dpfg[frame].PixelMeasuresSequence[0].PixelSpacing
    IPPMV = dpfg[frame].PlanePositionSequence[0].ImagePositionPatient
    IOPMV = dpfg[frame].PlaneOrientationSequence[0].ImageOrientationPatient

    PS = np.array([float(PSMV[0]), float(PSMV[1])])
    IPP = np.array([float(IPPMV[0]), float(IPPMV[1]), float(IPPMV[2])])
    IOP = np.array([float(IOPMV[0]), float(IOPMV[1]), float(IOPMV[2]), float(IOPMV[3]), float(IOPMV[4]), float(IOPMV[5]),])

    print(" ")
    print("Frame: " + str(frame+1) + ". Pixel Spacing (mm): " + str(PS))
    print("IPP: " + str(IPP))
    print("  IOP1: " + str(IOP[0:3]) )
    print("  IOP2: " + str(IOP[3:6]) )

    data = ds.pixel_array
    print("Array has shape: " + str(data.shape) + " and Fortran order is: " + str(data.flags.f_contiguous))


    img = data[frame] # bit like MATLAB data(frame,:,:) as frame is 1st dim here

    dshape = img.shape 
    nrow = dshape[0]
    ncol = dshape[1]

    L = np.zeros((nrow+1, ncol+1))
    P = np.zeros((nrow+1, ncol+1))
    H = np.zeros((nrow+1, ncol+1))


    if method=="vector addition":
        
        #  MATLAB:  TLV = IPP + -0.5*(IOP(1:3)*PS(2) + IOP(4:6)*PS(1))  ;
        TLV = IPP + -0.5*(IOP[0:3]*PS[1] + IOP[3:6]*PS[0])

        for irow in range(0,nrow+1):
            for icol in range(0, ncol+1):
                coordv = TLV + (irow)*IOP[3:6]*PS[0] + (icol)*IOP[0:3]*PS[1] ;
        
                L[irow, icol] = coordv[0]
                P[irow, icol] = coordv[1]
                H[irow, icol] = coordv[2]
                
    elif method=="matrix 2Dim":
                
        A = IOP[0:3] ; B = IOP[3:6] ;
        
        TM = np.array([  [ B[0], A[0], IPP[0] ], 
                         [ B[1], A[1], IPP[1] ],
                         [ B[2], A[2], IPP[2] ],
                         [ 0 ,   0,    1      ]  ] )
        
        S = np.array([ [ PS[0],    0,   0 ],
                       [ 0,    PS[1],   0 ],
                       [ 0,        0,   1 ] ] )
        
           
        for ir in range(0,nrow+1):
            for ic in range(0,ncol+1):
                v1 = ir-0.5  # image row coordinate for vertex
                v2 = ic-0.5  # image row coordinate for vertex 
                              # Note Python TLV will be at (-0.5 -0.5)
                
                # Multiply transform matrices with image coordinate
                # Python default has top left centre at (0,0).
                VLPH = np.matmul( np.matmul(TM, S), np.array([ [v1], [v2], [1]]))
                
                L[ir,ic] = VLPH[0] #  separate the components for surf plot
                P[ir,ic] = VLPH[1] 
                H[ir,ic] = VLPH[2] 
 
        
    else:
        print("Unsupported method supplied")
        
    
    scamap = plt.cm.ScalarMappable(cmap='gray')
    fcolors = scamap.to_rgba(img, alpha=0.5)
    ax.plot_surface(L, P, H, facecolors=fcolors, cmap='gray', linewidth=0, rcount=100, ccount=100)
    ax.set_xlabel('Left')
    ax.set_ylabel('Posterior')
    ax.set_zlabel('Head')
    
# end ddisp   -------------------



fig = plt.figure()
ax = fig.gca()

dslice_info(ax, fn_cor, fpath)
dslice_info(ax, fn_tra, fpath)
dslice_info(ax, fn_sag, fpath)
dslice_info(ax, fn_obl, fpath)

fig = plt.figure()
ax = fig.gca()

dslice(ax, fn_cor, fpath, 14, 1)
dslice(ax, fn_tra, fpath, 2, 2)
dslice(ax, fn_sag, fpath, 22, 3)
dslice(ax, fn_obl, fpath, 19, 4)

# For 3D plots, use method = "vector addition"  or method = "matrix 2Dim" 
fig = plt.figure()
ax = plt.axes(projection='3d')

method = "matrix 2Dim"
ddisp(ax, fn_cor, fpath, method, 14)
ddisp(ax, fn_tra, fpath, method,  2)
ddisp(ax, fn_sag, fpath, method, 22)
ddisp(ax, fn_obl, fpath, method, 19)


fig = plt.figure()
ax = plt.axes(projection='3d')

method = "vector addition"
ddisp(ax, fn_cor, fpath, method, 14)
ddisp(ax, fn_tra, fpath, method,  2)
ddisp(ax, fn_sag, fpath, method, 22)
ddisp(ax, fn_obl, fpath, method, 19)









