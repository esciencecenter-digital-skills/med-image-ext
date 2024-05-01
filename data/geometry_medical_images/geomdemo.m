function geomdemo
% GEOMDEMO Demonstrates geometry for DICOM files.
%  Accompanies "Geometry in Medical Imaging: DICOM and NIfTI Formats"
%  See also Python version  geomdemo.py
%
% Data available from Zenodo:
%  https://doi.org/10.5281/zenodo.6466491
%
% Run in the downloaded folder or edit the fpath variable 
% to point to folder with the DICOM files.
%
% Copyright 2022. David Atkinson, D.Atkinson@ucl.ac.uk
% University College London
%

% fpath is the full path to the folder containing the DICOM files IM_0002 etc
fpath = fullfile(pwd,'DICOM') ;

% File names
fn_cor = 'IM_0002' ; % Coronal
fn_tra = 'IM_0005' ; % Transverse
fn_sag = 'IM_0008' ; % Sagittal
fn_obl = 'IM_0011' ; % Oblique

% Check the slice information
figure('Name','Slice Info')
dslice_info(fullfile(fpath,fn_cor), fn_cor ) 
dslice_info(fullfile(fpath,fn_tra), fn_tra ) 
dslice_info(fullfile(fpath,fn_sag), fn_sag ) 
dslice_info(fullfile(fpath,fn_obl), fn_obl ) 

% Display specified slices in 3D 
method = 'vector addition' ;
figure('Name',method)
ddisp(fn_cor, fpath, method, 15)
ddisp(fn_tra, fpath, method,  3)
ddisp(fn_sag, fpath, method, 23)
ddisp(fn_obl, fpath, method, 20)

method = 'matrix 2Dim' ;
figure('Name',method)
ddisp(fn_cor, fpath, method, 15)
ddisp(fn_tra, fpath, method,  3)
ddisp(fn_sag, fpath, method, 23)
ddisp(fn_obl, fpath, method, 20)

figure('Name','Slices')
dslice(fn_cor, fpath, 15, 1)
dslice(fn_tra, fpath, 3, 2)
dslice(fn_sag, fpath, 23, 3)
dslice(fn_obl, fpath, 20,4)

end


function dslice(fn, fpath, frame, splot)
% Displays one slice on a 2x2 plot
ffn = fullfile(fpath, fn) ;
data = dicomread(ffn) ;
img = data(:,:,frame) ;
subplot(2,2,splot)
imshow(img,[])
title(fn,'Interpreter','none')
end


function dslice_info(ffn, label)
% DSLICE_INFO Display Dicom slice information 
%  Plots the slice position along the slice normal.
%  Checks for parallel slices.
%  Checks first two slices are in correct order.
%  Commented code provided for sorting slices.
% 
% dslice_info(ffn, label) 
%  ffn    full file name including any path needed
%  label  user-supplied label, typically the file name 
%

[~, spatialInfo, ~] = dicomreadVolume(ffn) ;
IOPs = spatialInfo.PatientOrientations ;  % [2 x 3 x nframe]
IPPs = spatialInfo.PatientPositions ;     % [nframe x 3]

nframe = size(IPPs,1) ;

% Use first frame's ImageOrientationPatient as reference and 
% compute the slice normal C

C = cross(IOPs(1,:,1), IOPs(2,:,1) ) ;  % slice normal

proj = zeros(1,nframe) ; % allocate space for projections

for iframe = 1: nframe
    % Projection of frame's ImagePositionPatient along slice normal
    proj(1,iframe) = dot(C, IPPs(iframe,:) );
    
    % Check ImagePositionPatient is unchanged (i.e. slices are parallel)
    if norm(IOPs(:,:,iframe) - IOPs(:,:,1)) > 0.001
        warning('Slices may not be parallel')
    end
end

plot(proj,'LineWidth',3,'DisplayName',label)

fs = 14 ; % FontSize
legend('interpreter','none', 'FontSize', fs)
title('Projections of IPP onto slice normal', 'FontSize', fs)
xlabel('Frame number', 'FontSize', fs), ylabel('Projection (mm)', 'FontSize', fs)
hold on, grid on

disp(['Distance between slice 1 and 2 centres is: ', ...
       num2str(proj(2) - proj(1)), ' mm for data: ',label])

if proj(2) < proj(1)
    warning(['First and 2nd slices are not in ascending order for: ',label])
end


% To be more complete, sort proj into ascending order and use the 
% sort indices to re-order the data and spatial information. 
% Also check that slice separations are equal.
%  
%  [sortedproj, idx] = sort(proj,2) ;
%  Vout = V(:,:,:,idx) ;
%  IPPsout = IPPs(idx,:) ;
%  IOPsout = IOPs(:,:,idx) ;
%
%  if norm(diff(diff(sortedproj))) > 0.01 
%    warning(['Slice separations are not equal'])
%  end
%
%
% dicomreadVolume documentation implies it does do this sorting, but it
% appears not to be the case for a multi-frame DICOM.

end

% - - - - - - - - - - - - - 

function ddisp(fn, fpath, method, frame)
% DDISP Dicom display on 3D plot
%  Displays slices on 3D plot.
%
%  ddisp(fn, fpath, method, frame)
%  ddisp(fn, fpath, method)  defaults to middle slice for volume plot
%
% fn is filename, fpath is path to this file
% method can be 'matrix 2Dim' or 'vector addition' 
%

ffn = fullfile(fpath, fn) ;

dinfo = dicominfo(ffn) ; 

% For more compact reading, use dicomreadVolume - see example in dslice_info
% above.
if strcmp(dinfo.SOPClassUID, '1.2.840.10008.5.1.4.1.1.4.1' )
    % enhanced MR (multi-frame)
    if nargin < 3
        frame = ceil(dinfo.NumberOfFrames/2) ;
    end
    fieldslice = ['Item_',num2str(frame)] ;
    dinfoslice = dinfo.PerFrameFunctionalGroupsSequence.(fieldslice) ;
    IOP = dinfoslice.PlaneOrientationSequence.Item_1.ImageOrientationPatient ;
    IPP = dinfoslice.PlanePositionSequence.Item_1.ImagePositionPatient ;
    PS  = dinfoslice.PixelMeasuresSequence.Item_1.PixelSpacing ;
    
    disp(['Frame: ',num2str(frame),', IPP: ',num2str(IPP')])
    disp(['  IOP: ',num2str(IOP')])
else
    % Code here is just for completeness - not used by this example
    frame = 1 ;
    IOP = dinfo.ImageOrientationPatient ;
    IPP = dinfo.ImagePositionPatient ;
    PS =  dinfo.PixelSpacing ;
end

dat = dicomread(dinfo) ; % read image data
vdat = squeeze(double(dat)) ;
%nframe = size(vdat,3) ;

img = mat2gray(vdat(:,:,frame)) ; % 2D image

[nrow, ncol] = size(img) ;

% Allocate matrices L, P and H to contain image vertex coordinates
L = zeros(nrow+1, ncol+1) ;
P = zeros(nrow+1, ncol+1) ;
H = zeros(nrow+1, ncol+1) ;

% 3D coord of top left vertex.
%
%  TLV
%  +-----+-----+-----+   IOP(1:3) ->
%  |     |     |     |   
%  |  x  |     |     |   |
%  |     |     |     |   v IOP(4:6)
%  +-----+-----+-----+
%  |     |     |     |   PS = [height width]
%  |     |     |     |
%
% x is centre of top left coord, IPP in the patient coordinate system.
% To compute TLV (top left vertex), go back half a pixel and up half a pixel.
% Scale for the pixel spacing in mm (PS).

switch method
    case 'vector addition'  % Vector Addition method
        
        TLV = IPP + -0.5*(IOP(1:3)*PS(2) + IOP(4:6)*PS(1))  ;
        
        % loop over each vertex
        for ir = 1:nrow+1
            for ic = 1:ncol+1
                % vertex coordinate, coordv, is just the coord of the
                % top left vertex, TLV, plus a vector down the rows
                % plus a vector along the columns.
                coordv = TLV + (ir-1)*IOP(4:6)*PS(1) + (ic-1)*IOP(1:3)*PS(2) ;
                
                L(ir,ic) = coordv(1) ; % separate the components for surf plot
                P(ir,ic) = coordv(2) ;
                H(ir,ic) = coordv(3) ;
            end
        end
        
    case 'matrix 2Dim'
        % 2D image coordinate and transformation matrix to calcualte vertex
        % coordinates
        
        A = IOP(1:3) ; B = IOP(4:6) ;  % A is IOP1 in report and B is IOP2
        
        TM = [ B(1) A(1) IPP(1)  ; % note order
               B(2) A(2) IPP(2)  ;
               B(3) A(3) IPP(3)  ;
               0    0    1     ] ;
        
        S = [PS(1) 0     0 ;
            0      PS(2) 0 ;
            0      0     1 ] ;
        
        Timtlc = [ 1 0 -1 ;    % Offset for top left centre image coordinate 
                   0 1 -1 ;    % This is (1,1) in MATLAB default system
                   0 0  1 ] ;
           
        for ir = 1:nrow+1
            for ic = 1:ncol+1
                v1 = ir-0.5 ; % image row coordinate for vertex
                v2 = ic-0.5 ; % image row coordinate for vertex 
                              % Note TLV will be at (0.5 0.5)
                
                % Multiply transform matrices with image coordinate
                %   (written as a homogeneous column vector)
                VLPH = TM * S * Timtlc * [v1; v2; 1];
                
                L(ir,ic) = VLPH(1) ; % separate the components for surf plot
                P(ir,ic) = VLPH(2) ;
                H(ir,ic) = VLPH(3) ;
            end
        end
        
    otherwise
        error(['Unknown method: ', method])
end

% Surf plot in 3D for 'vector addition' or 'matrix 2Dim' methods

surf(L,P,H,img,'EdgeColor','None', 'FaceAlpha',0.5)

hold on, grid on
axis equal
axis vis3d
view(21,-50)
xlabel('Left') ; ylabel('Posterior') ; zlabel('Head')
colormap gray

end
