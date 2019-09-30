function defField = calcDefFieldSliding(cpg1_nii, cpg2_nii, dist12_nii, target_nii, rxs, rys, rzs )
%calcDefFieldSliding function to calculate a deformation field from two 
% B-spline control point grids and a signed distance function according to
% the sliding b-Spline registration approach
%
%INPUTS: cpg1_nii:   a nifti structure of a control-point grid where the 
%                    image data contains a 4D array describing the B-spline 
%                    control point displacements for region 1 i.e. 
%                    transformed and added signed distance function < 0
%        cpg2_nii:   same as CPG1_nii, but describing transformation of 
%                    region 2, i.e. transformed and added signed distance 
%                    function >= 0 
%        dist12_nii: signed distance functions (warped and added in a pre-
%                    processing step) in the target image space identifying 
%                    sliding regions (reg1 < 0, reg2 >= 0)
%        target_nii: nifti image strucure holding the targt image (required 
%                    to transform the distance map into the target image 
%                    space)
%        rxs: the x coords of the reference image (i.e. the points where
%        the deformation field will be calculated)
%        rys: the y coords of the reference image
%        rzs: the z coords of the reference image
%OUTPUTS: defField: a 4D array containing the deformation field
% NOTE: dist12

%check for optional inputs
if ~exist('rxs','var') 
    [rxs,rys,rzs] = coords_from_nii(target_nii);
end


% Calculate two deformation vector fields, as if both control point grids
% defined a separate transformation
[cxs,cys,czs] = coords_from_nii(cpg1_nii);
defField1 = calcDefField( permute(cpg1_nii.img, [1,2,3,5,4]), ...
                          cxs,cys,czs,rxs,rys,rzs);
defField2 = calcDefField( permute(cpg2_nii.img, [1,2,3,5,4]), ...
                          cxs,cys,czs,rxs,rys,rzs);
                        
% the distance function coordinates
[dxs, dys, dzs] = coords_from_nii( dist12_nii );
[dxsG,dysG,dzsG] = ndgrid( dxs, dys, dzs );
[rxsG,rysG,rzsG] = ndgrid( rxs, rys, rzs );

% we need to determine, which points of the reference coordinates belong to
% region 1 and which belong to region 2
% --> sample the added distance map (dist12) at the specified locations
if cpg1_nii.hdr.dime.dim(6)==2
  distMapSamples = interpn( dxsG, dysG, double( dist12_nii.img ), rxsG, rysG);
  distMapSamples = cat(4, distMapSamples, distMapSamples);
else
  distMapSamples = interpn( dxsG, dysG, dzsG, double( dist12_nii.img ), rxsG, rysG, rzsG);
  distMapSamples = cat(4, distMapSamples, distMapSamples, distMapSamples);
end
% Compose the deformation field
% Region 1 <  0
% Region 2 >= 0
% NaN everywehere else
defField = nan(size(defField1));
defField(distMapSamples<0)  = defField1(distMapSamples <  0);
defField(distMapSamples>=0) = defField2(distMapSamples >= 0);

end