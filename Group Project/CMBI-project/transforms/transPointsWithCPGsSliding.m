function points_trans = transPointsWithCPGsSliding(cpg1_nii, cpg2_nii, dist_nii, points_orig, target_nii, cpg_is_disp, reorientate)
%function to transform one or more points using the provided control point
%grids  which defines a b-spline transformation obtained using sliding
%b-spline registration 
%
%INPUTS:
%   cpg1_nii - the control point grid as a nifti structure for region 1
%   cpg2_nii - the control point grid as a nifti structure for region 2
%   dist_nii - the signed distance function in the source image space
%   points_orig - an n x 3 vector (or n x 2 if 2D) containing the point
%   coordinates before applying the transformtion
%   cpg_is_disp - if true the cpg contains a displacement field, if false
%       the cpg contains a deformation field. [false]
%   reorientate - correct inversion of first and second dimension. [true]
%OUTPUTS:
%   points_trans - an n x 3 vector containing the point coordinates after
%   applying the transformtion
%
%Note: deformation field = target coords + displacement field
%Note: For the CMBI-project the input parameter reorientate defaults to
%      true
%check for optional inputs
if ~exist('cpg_is_disp','var') || isempty(cpg_is_disp)
    cpg_is_disp = false;
end

if ~exist('reorientate','var') || isempty(reorientate)
    reorientate = true;
end

if reorientate
    points_orig(:,1:2) = -points_orig(:,1:2);
end

%get coords for cpg
[cxs,cys,czs] = coords_from_nii(cpg1_nii);

%loop over points and apply transformation
points_trans = zeros(size(points_orig));
if size(points_orig,2) == 2
    points_orig = [points_orig zeros(size(points_orig,1),1)];
end

% Pre-calculate the warped 
dist12_nii = preCalculateDeformedSignedDistForSlidingBSplineDefField(cpg1_nii, cpg2_nii, dist_nii, target_nii, cpg_is_disp );
parfor pt = 1:size(points_orig,1)
    %get deformation/displacement for this point
    points_trans(pt,:) = squeeze( calcDefFieldSliding(cpg1_nii, cpg2_nii, dist12_nii, target_nii, points_orig(pt,1),points_orig(pt,2),points_orig(pt,3) ) );
end

%if cpg contains displacements need to add original point coords
if cpg_is_disp
    points_trans = points_trans + points_orig;
end
if reorientate
    points_trans(:,1:2) = -points_trans(:,1:2);
end

end