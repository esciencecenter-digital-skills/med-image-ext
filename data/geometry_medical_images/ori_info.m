function oinfo = ori_info(vec, varargin)
% ORI_INFO Get orientation information.
% Determines the Radiological name of the plane with the closest normal 
% to the input. 
% Input is either a 3-element vector, or a 6-element variable (corrresponding
% to ImageOrientationPatient) where the cross-product of the two 3-element
% vectors is used. 
% Also returns the east, west, north and south strings if ImageOrientationPatient
% was input.
% Code does not check if orientation meets radiological viewing conventions.
%
%  oinfo = ori_info(vec_lph)  3 element vector in LPH system (slice normal)
%  oinfo = ori_info(IOP)      6 element IOP array
%  oinfo = ori_info(..., Name, Value,...)
%
% Name Value pairs:
%  'label_tolerance' {0.01}  dot product tolerance for west_str etc
%  'geom' geom structure     checks through plane
%
%  oinfo is a structure with fields:
%    plane_str           'TRA', 'SAG', or 'COR' 
%    closest_normal_str  'RL', 'LR', 'AP', 'PA', 'FH' or 'HF' 
%                         based on cross product of IOP
%    If geom was input:
%        positive_z_str        needs geom input
%       isRH      true, false  need geom input
%
%    If IOP was input:
%         west_str, east_str, north_str, south_str  Labels for sides of images.
%         Examples: 'L'   if image plane is aligned to within a tolerance
%                   'LA'  if single angulation
%                   'LAH' if double oblique
%
% Copyright, 2020-2022, David Atkinson, University College London
% D.Atkinson@ucl.ac.uk
%
% Copied from tools for geom paper.
%
% See also NDCOMP

label_tolerance = 0.01 ;
geom = [] ;
isRH = [];

for ipv = 1:2:length(varargin)
    val = varargin{ipv+1} ;
    switch varargin{ipv}
        case 'label_tolerance'
            label_tolerance = val ;
        case 'geom'
            geom = val ;
        otherwise
            error(['Unknown input: ',varargin{ipv}]) ;
    end    
end
                
if length(vec) == 6
    sdc = cross(vec(1:3), vec(4:6)) ;
elseif length(vec) == 3
    sdc = vec ;
else
    error('First arg must be normal vector or IOP vectors')
end

if abs(norm(sdc)-1) > 0.01
    warning(['Expecting vector norm of 1, found: ', num2str(norm(sdc)),' Normalising'])
    sdc = sdc ./ norm(sdc) ;
end

        
sagnorm = cross([0 1 0],[0 0 -1]) ; % LPH coordinate
cornorm = cross([1 0 0],[0 0 -1]) ;
tranorm = cross([1 0 0],[0 1  0]) ;

sagproj = dot(sagnorm, sdc) ; 
corproj = dot(cornorm, sdc) ;
traproj = dot(tranorm, sdc) ;

projs = [traproj sagproj corproj ] ;
oristr  = {'TRA', 'SAG', 'COR'} ; 
    
[~, iproj] = max(abs(projs));

vds = [ 1 0 0 ; -1 0 0 ; 0 1 0 ; 0 -1 0 ; 0 0 1 ; 0 0 -1] ;
dstr = { 'RL' , 'LR'   , 'AP' , 'PA' , 'FH' , 'HF' } ;

dps = dot(repmat(sdc(:)',[6 1]), vds, 2) ;
[~, imdp] = max(dps) ;

% disp(['Vec is ',dstr{imdp},'. Plane norm with closest orientation is: ',...
%        oristr{iproj}, ' having dot product: ', num2str(projs(iproj))])

% This if-statement is not relevant to the Geometry paper
if ~isempty(geom)
    if length(geom)> 1
        tsd = geom(2).IPP - geom(1).IPP ;
        tsd = tsd ./norm(tsd) ;
        
        tsddp = dot(repmat(tsd(:)',[6 1]), vds, 2) ;
        [~, imtsddp] = max(tsddp) ;
        
        oinfo.positive_z_str = dstr{imtsddp} ;
        
        if abs(dot(sdc, tsd) - 1 ) < 0.01
            isRH = true ;
        end
        if abs(dot(sdc, tsd) + 1 ) < 0.01
            isRH = false ;
        end
          
        oinfo.isRH = isRH ;
    end
end


oinfo.plane_str = oristr{iproj} ;
oinfo.closest_normal_str  = dstr{imdp} ;
if length(vec) == 6 % it was an IOP
    iop = vec ;
    oinfo.west_str = determine_label(-iop(1:3), label_tolerance) ;
    oinfo.east_str = determine_label(iop(1:3), label_tolerance) ;
    oinfo.north_str = determine_label(-iop(4:6), label_tolerance) ;
    oinfo.south_str = determine_label(iop(4:6), label_tolerance) ;
end

end


function label_str = determine_label(vec_lph, label_tolerance)
% DETERMINE_LABEL Determines the label for an LPH direction
%  Used to label the edges of an image
%
% label_str = determine_label(vec_lph, label_tolerance)
%
%  vec_lph  - vector in LPH
% 
%  label_str label for direction of vec_lph
% 

% Find predominant direction
ax_vec = [ 1 0 0 ; -1 0 0 ; 0 1 0 ; 0 -1 0 ; 0 0 1 ; 0 0 -1] ;
ax_str = { 'L' , 'R'   ,  'P' ,    'A' ,    'H' ,   'F' } ;

dotax = dot(repmat(vec_lph(:)',[6 1]), ax_vec, 2) ;
[sorted_dotax, ind_sorted] = sort(dotax,'descend') ;

label_str(1) = ax_str{ind_sorted(1)} ;

if sorted_dotax(2) > label_tolerance
    label_str(2) = ax_str{ind_sorted(2)} ;
    if sorted_dotax(3) > label_tolerance
        label_str(3) = ax_str{ind_sorted(3)} ;
    end
end


end



