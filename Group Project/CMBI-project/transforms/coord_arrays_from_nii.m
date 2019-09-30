function [X, Y, Z] = coord_arrays_from_nii(nii)
%function to get the x, y, and z voxel coordinates from a nifti structure
%and return them as 3D arrays
%
%currently only works for sform

if nii.hdr.hist.sform_code ~= 1
    error('sform must be set to calculate coordinate arrays');
else
    
    %3D arrays of voxel indices
    [Xv, Yv, Zv] = ndgrid(0:nii.hdr.dime.dim(2)-1,0:nii.hdr.dime.dim(3)-1,0:nii.hdr.dime.dim(4)-1);
    
    %calc arrays with mm coords using sform matrix
    X = Xv * nii.hdr.hist.srow_x(1) + Yv * nii.hdr.hist.srow_x(2) + Zv * nii.hdr.hist.srow_x(3) + nii.hdr.hist.srow_x(4);
    Y = Xv * nii.hdr.hist.srow_y(1) + Yv * nii.hdr.hist.srow_y(2) + Zv * nii.hdr.hist.srow_y(3) + nii.hdr.hist.srow_y(4);
    Z = Xv * nii.hdr.hist.srow_z(1) + Yv * nii.hdr.hist.srow_z(2) + Zv * nii.hdr.hist.srow_z(3) + nii.hdr.hist.srow_z(4);
    
end