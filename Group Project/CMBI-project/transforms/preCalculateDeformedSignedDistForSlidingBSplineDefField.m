function [dist12_nii, dist1_nii, dist2_nii] = preCalculateDeformedSignedDistForSlidingBSplineDefField(cpg1_nii, cpg2_nii, dist_nii, target_nii, cpg_is_disp )
%
% Function to speed transformation of points, when the deformation field is
% repeatedly calculated
%

  %check for optional inputs
  if ~exist('cpg_is_disp','var') || isempty(cpg_is_disp)
      cpg_is_disp = false;
  end

  % This was placed outside the function to calculate the defFieldSpliding
  dist1_nii = deformNiiWithCPG(cpg1_nii, dist_nii, target_nii, cpg_is_disp);
  dist2_nii = deformNiiWithCPG(cpg2_nii, dist_nii, target_nii, cpg_is_disp);

  % Add the distance maps
  dist12_nii = dist1_nii;
  dist12_nii.img = dist1_nii.img + dist2_nii.img;
end