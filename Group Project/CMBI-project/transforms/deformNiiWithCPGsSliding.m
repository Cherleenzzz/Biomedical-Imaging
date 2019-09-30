function [def_vol_nii, def_field_nii, dis_field_nii] = deformNiiWithCPGsSliding(cpg1_nii, cpg2_nii, dist_nii, source_nii, target_nii, cpg_is_disp)
% function to deform a source image into the space of a target image using
% the provided control point grids and signed distance image representing a 
% sliding b-spline transformation. Note that the interface differs from the
% corresponding non-sliding matlab function.
%
%INPUTS:
%   cpg1_nii - the control point grid as a nifti structure for region 1
%   cpg2_nii - the control point grid as a nifti structure for region 2
%   dist_nii - the signed distance function in the source image space
%              defining which parts belong to region 1 and which to 
%              region 2
%   source_nii - the source image as a nifti structure
%   target_nii - the target image as a nifti structure
%   cpg_is_disp - if true the cpg contains a displacement field, if false
%       the cpg contains a deformation field. [false]
%OUTPUTS:
%   def_vol_nii - the deformed image as a nifti structure
%   def_field_nii - the deformation field as a nifti structure
%   dis_field_nii - the displacement field as a nifti structure
%
%Note: deformation field = target coords + displacement field

%images and cpg can be 2D or 3D
%number of components (2 or 3) to vectors in CPG used to determine number
%of dimensions
is3D = true;
if cpg1_nii.hdr.dime.dim(6) == 2
    is3D = false;
elseif cpg1_nii.hdr.dime.dim(6) ~= 3
    msg = 'CPG must contain 2D or 3D vectors\n';
    msg = [msg 'these should be stored on along the 5th dimesnion.'];
    error(msg);
end

%check for optional inputs
if ~exist('cpg_is_disp','var') || isempty(cpg_is_disp)
    cpg_is_disp = false;
end

% get coords for images
[sxs,sys,szs] = coords_from_nii(source_nii);
[txs,tys,tzs] = coords_from_nii(target_nii);

% form matrices from image coords using ndgrid
[sX,sY,sZ] = ndgrid(sxs,sys,szs);
[tX,tY,tZ] = ndgrid(txs,tys,tzs);

% calculate def field from cpg
% pre-caluculate the signed distnance map (warped and added)
dist12_nii = preCalculateDeformedSignedDistForSlidingBSplineDefField(cpg1_nii, cpg2_nii, dist_nii, target_nii, cpg_is_disp );
def_field  = calcDefFieldSliding(cpg1_nii, cpg2_nii, dist12_nii, target_nii);

if cpg_is_disp
    if is3D
        def_field = def_field + cat(4,tX,tY,tZ);
    else
        def_field = def_field + cat(4,tX,tY);
    end
end
    
%deform source image using def_field
def_vol_nii = target_nii;

%check for 2D or 3D
if is3D
    def_vol_nii.img = interpn(sX,sY,sZ,double(source_nii.img),def_field(:,:,:,1),def_field(:,:,:,2),def_field(:,:,:,3),'spline', nan);
else
    def_vol_nii.img = interpn(sX,sY,double(source_nii.img),def_field(:,:,:,1),def_field(:,:,:,2));
end

if nargout > 1
    def_field_nii = target_nii;
    def_field_nii.hdr.dime.dim(1) = 5;
    def_field_nii.hdr.dime.intent_p1 = 0; % Follow nifty-reg conventions
    if is3D
        def_field_nii.hdr.dime.dim(6) = 3;
    else
        def_field_nii.hdr.dime.dim(6) = 2;
    end
    def_field_nii.hdr.dime.datatype = cpg1_nii.hdr.dime.datatype;
    def_field_nii.hdr.dime.bitpix = cpg1_nii.hdr.dime.bitpix;
    df_size = size(def_field);
    def_field_nii.img = reshape(def_field,df_size(1),df_size(2),df_size(3),1,df_size(4));
    if nargout > 2
        dis_field_nii = def_field_nii;
        dis_field_nii.hdr.dime.intent_p1 = 1; % Follow nifty-reg conventions
        if is3D
            dis_field_nii.img = def_field - cat(4,tX,tY,tZ);
        else
            dis_field_nii.img = def_field - cat(4,tX,tY);
        end            
    end 
end
