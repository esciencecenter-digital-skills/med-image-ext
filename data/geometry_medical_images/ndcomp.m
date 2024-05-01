% NDCOMP Compares transforms from NIfTI and DICOM 
%
% Accompanies paper "Geometry in Medical Imaging: DICOM and NIfTI formats"
%
% Run from folder that contains DICOM and NIFTI folders.
%
% Copyright 2022. David Atkinson, University College London
% D.Atkinson@ucl.ac.uk
%
% See also ORI_INFO GEOMDEMO


compt('DICOM/IM_0002', 'NIFTI/OBJECT_phantom_T2W_TSE_Cor_14_1.nii', 15);
compt('DICOM/IM_0005', 'NIFTI/OBJECT_phantom_T2W_TSE_Tra_17_1.nii', 3);
compt('DICOM/IM_0008', 'NIFTI/OBJECT_phantom_T2W_TSE_Sag_18_1.nii', 23 );
compt('DICOM/IM_0011', 'NIFTI/OBJECT_phantom_T2W_TSE_obl_19_1.nii', 20) ;


compt('DICOM/IM_0002', 'NIFTI/OBJECT_phantom_T2W_TSE_Cor_FSL_14_1.nii', 15);
compt('DICOM/IM_0005', 'NIFTI/OBJECT_phantom_T2W_TSE_Tra_FSL_17_1.nii', 3);
compt('DICOM/IM_0008', 'NIFTI/OBJECT_phantom_T2W_TSE_Sag_FSL_18_1.nii', 23);
compt('DICOM/IM_0011', 'NIFTI/OBJECT_phantom_T2W_TSE_obl_FSL_19_1.nii', 20) ;


function compt(dfile, nfile, dslice)
% COMPT Compare DICOM and NIfTI transforms and display with anatomic
% direction
%
% DICOM transform is built using IPP from first slice (Item_1)
%
% Copyright 2020-2022. David Atkinson, University College London
% D.Atkinson%ucl.ac.uk
%

di = dicominfo(dfile);
IPP = di.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient ;
IOP = di.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient;

% DICOM IOP is [ unit_vec_along_row  unit_vec_down_column] 

N = cross(IOP(1:3), IOP(4:6)) ;  % slice normal (C in paper)

% patient3D_coord = T * img_coord
%  where patient2D_coord and img_coord are column vectors in homogeneous coords
%
% Row 1 of T corresponds to L (DICOM) or R (NIfTI).
% Row 2 of T corresponds to P (DICOM) or A (NIfTI).
% Row 3 of T corresponds to H (DICOM) or H (NIfTI).
%  For img_coord(i1, i2, i3):
% Col 1 of T is a unit vector in the directon (LPH/RAH) of increasing i1.
% Col 2 of T is a unit vector in the directon (LPH/RAH) of increasing i2.
% Col 3 of T is a unit vector in the directon (LPH/RAH) of increasing i3.

% In MATLAB default image coordinates i1 is row-number and increasing i1
% goes down a column  - this is IOP(4:6)
% Increasing i2 goes left to right along a row - IOP(1:3)

% Find from DICOM, transform expected for NIfTI
% Negate the first two rows as NIfTI is RAH wheras DICOM is LPH
T_RAH = [-IOP(4) -IOP(1) -N(1)  -IPP(1) ;
         -IOP(5) -IOP(2) -N(2)  -IPP(2) ;
          IOP(6)  IOP(3 ) N(3)   IPP(3) ;
           0       0        0       1    ] ;

ni = niftiinfo(nfile) ;
nit = ni.Transform.T' ; % NIfTI transform

% Account for scaling that is included in nit
nit(:,1) = nit(:,1)/ni.PixelDimensions(1) ;
nit(:,2) = nit(:,2)/ni.PixelDimensions(2) ;
nit(:,3) = nit(:,3)/ni.PixelDimensions(3) ;

% Determine DICOM-style IOP from NIfTI transform (without scaling)
% negative signs for rows 1&2 are to convert from NIfTI RAH to DICOM LPH.
IOPfromNIT = [-nit(1,2) -nit(2,2) nit(3,2) -nit(1,1) -nit(2,1) nit(3,1)] ;

% Determine if transforms are equal. Usually they are not because the
% images are rotated with respect to each other. However, on examining the
% displayed images, you can check that the directions correctly correspond.
diffT = nit - T_RAH ;

if max(abs(diffT(:))) > 1e-5
    disp(['Transforms differ between ',dfile,' and ',nfile])
else
    disp(['Tranforms are the same from ',dfile,' and ',nfile])
end

% Show dslice from each and label the patient direction
figure('Name',[dfile,' and ',nfile])
dv = dicomread(di) ;
niv = niftiread(ni) ;

t = tiledlayout(1,2);
fs = 14 ;

nexttile
imshow(dv(:,:,1,dslice),[])
title(dfile,'Interpreter','none', 'FontSize',fs)
oinfo = ori_info(IOP, 'label_tolerance',0.1) ;
xlabel(oinfo.south_str, 'FontSize',fs)
ylabel(oinfo.west_str, 'FontSize',fs)


nexttile
imshow(niv(:,:,dslice),[])
title(nfile,'Interpreter','none', 'FontSize',fs)
oinfo = ori_info(IOPfromNIT, 'label_tolerance',0.1) ;
xlabel(oinfo.south_str, 'FontSize',fs)
ylabel(oinfo.west_str, 'FontSize',fs)

end
 
 
 
