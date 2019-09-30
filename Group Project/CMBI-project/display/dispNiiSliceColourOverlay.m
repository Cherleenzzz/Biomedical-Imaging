function h = dispNiiSliceColourOverlay(nii1,nii2,dim,slice,win1,win2,cmap1,cmap2,switch_dims,reorientate,use_vox_coords,show_axes)
%h = dispNiiSliceColourOverlay(nii1,nii2,dim,slice,win1,win2,cmap1,cmap2,switch_dims,reorientate,use_vox_coords)
%function to display a slice from two nifti volumes using a colour overlay
%
%INPUTS: nii1: nifti structure containing the 1st image and header info
%       nii2: nifti structure containing the 2nd image and header info
%       dim: the dimension that the slice is taken from ('x', 'y', or 'z')
%       slice: the slice number to display
%       win1: a two element vector with the intensity window to use when
%               displaying the 1st image [full intensity range]
%       win2: a two element vector with the intensity window to use when
%               displaying the 2nd image [full intensity range]
%       cmap1: the color map to use when displaying the 1st image [red]
%       cmap2: the color map to use when displaying the 2nd image [cyan]
%       switch_dims: boolean specifying whether x and y dimensions should
%               be switched before displaying [false]
%       reorientate: if true use -ve of x and y coords - this displays 
%               images same as in NiftyView (default = true if sform set,
%               false otherwise)
%       use_vox_coords: if true display using voxel coordinates rather than
%               world coordinates [false]
%
%OUTPUTS:   h: handle to the image object
%
%assumes the two nifti volumes have the same voxel coordinates

%check for optional inputs
if ~exist('win2','var') || isempty(win2)
    if ~exist('win1','var') || isempty(win1)
        win1 = [min(nii1.img(:)) max(nii1.img(:))];
        win2 = [min(nii2.img(:)) max(nii2.img(:))];
    else
        win2 = win1;
    end
end
if ~exist('cmap1','var') || isempty(cmap1)
    cmap1 = gray;
    cmap1(:,2:3) = 0;
end
if ~exist('cmap2','var') || isempty(cmap2)
    cmap2 = gray;
    cmap2(:,1) = 0;
end
if ~exist('switch_dims','var') || isempty(switch_dims)
    switch_dims = false;
end
if ~exist('reorientate','var') || isempty(reorientate)
    if nii1.hdr.hist.sform_code == 0
        reorientate = false;
    else
        reorientate = true;
    end
end
if ~exist('use_vox_coords','var') || isempty(use_vox_coords)
    use_vox_coords = false;
end
if ~exist('show_axes','var')|| isempty(show_axes)
    show_axes=true;
end

%get vox coords
if use_vox_coords
    if nii1.hdr.dime.dim(1:4) ~= nii2.hdr.dime.dim(1:4)
        error('images must have the same dimensions');
    end
    xs = 0:nii1.hdr.dime.dim(2)-1;
    ys = 0:nii1.hdr.dime.dim(3)-1;
    zs = 0:nii1.hdr.dime.dim(4)-1;
else
    [xs, ys, zs] = coords_from_nii(nii1);
    [xs2, ys2, zs2] = coords_from_nii(nii2);
    if any(xs2 ~= xs) || any(ys2 ~= ys) || any(zs2 ~= zs)
        error('images must have the same voxel coordinates');
    end
    if reorientate
        xs = -xs;
        ys = -ys;
    end
    %check for single slices
    if length(xs) == 1
        xs = xs + [-1 1]*nii1.hdr.dime.pixdim(2)/6;
    end
    if length(ys) == 1
        ys = ys + [-1 1]*nii1.hdr.dime.pixdim(3)/6;
    end
    if length(zs) == 1
        zs = zs + [-1 1]*nii1.hdr.dime.pixdim(4)/6;
    end
end


%select slice and set coords for slice
switch lower(dim)
    case 'x'
        im1 = squeeze(nii1.img(slice,:,:))';
        im2 = squeeze(nii2.img(slice,:,:))';
        xs_im = ys;
        ys_im = zs;
    case 'y'
        im1 = squeeze(nii1.img(:,slice,:))';
        im2 = squeeze(nii2.img(:,slice,:))';
        xs_im = xs;
        ys_im = zs;
    case 'z'
        im1 = nii1.img(:,:,slice)';
        im2 = nii2.img(:,:,slice)';
        xs_im = xs;
        ys_im = ys;
    otherwise
        error('axis must be x, y, or z');
end

if switch_dims
    tmp = xs_im;
    xs_im = ys_im;
    ys_im = tmp;
    im1 = im1';
    im2 = im2';
end


%convert im1(im2) into indexed image into cmap1(cmap2) using win1(win2)
im1 = gray2ind(mat2gray(im1,double(win1)),size(cmap1,1));
im2 = gray2ind(mat2gray(im2,double(win2)),size(cmap2,1));

%make true colour image using cmap1 and cmap2 and display
im = reshape(cmap1(im1+1,:)+cmap2(im2+1,:),[size(im1) 3]);
%calculate max combined colour and renormalise if needed
cmap_max = max(cmap1(end,:) + cmap2(end,:));
if cmap_max>1
    im = im/cmap_max;%rescale intensities
end
image(xs_im,ys_im,im);

%set axis directions
if ~isequal(dim,'z')
    set(gca,'ydir','normal');
end

%scale axis if using vox coords
axis image
if use_vox_coords
    switch lower(dim)
        case 'x'
            daspect([1./nii1.hdr.dime.pixdim(3:4) 1]);
        case 'y'
            daspect([1./nii1.hdr.dime.pixdim([2 4]) 1]);
        case 'z'
            daspect([1./nii1.hdr.dime.pixdim(2:3) 1]);
    end
end

if ~show_axes
    ax=gca();
    ax.Visible=0;
end
