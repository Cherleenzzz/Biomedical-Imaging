function h = dispNiiSlice(nii,dim,slice,clim,cmap,switch_dims,reorientate,use_vox_coords,show_axes)
%h = dispNiiSlice(nii,dim,slice,clim,cmap,switch_dims)
%function to display a slice of a 3D volume stored as a nifti structure
%
%INPUTS: nii: nifti structure containing image and header info
%       dim: the dimension that the slice is taken from ('x', 'y', or 'z')
%       slice: the slice number to display
%       clim: a two element vector with the min and max intensities to
%               display (default [min(V(:)) max(V(:))]
%       cmap: the color map to use when displaying the image (default
%               gray(256))
%       switch_dims: boolean specifying whether x and y dimensions should
%               be switched before displaying
%       reorientate: if true use -ve of x and y coords - this displays 
%               images same as in NiftyView (default = true if sform set,
%               otherwise false)
%       use_vox_coords: if true display using voxel coordinates rather than
%               world coordinates (default = false)
%
%OUTPUTS:   h: handle to the image object

V = nii.img;
%check for inputs
if ~exist('clim','var') || isempty(clim)
    clim = [min(V(:)) max(V(:))];
end
if ~exist('cmap','var') || isempty(cmap)
    cmap = gray(256);
end
if ~exist('switch_dims','var') || isempty(switch_dims)
    switch_dims = false;
end
if ~exist('reorientate','var') || isempty(reorientate)
    
    if nii.hdr.hist.sform_code > 0
        reorientate = true;
    else
        reorientate = false;
    end
end
if ~exist('use_vox_coords','var') || isempty(use_vox_coords)
    use_vox_coords = false;
end
if ~exist('show_axes','var') || isempty(show_axes)
    show_axes=true;
end

%get vox coords
if use_vox_coords
    xs = 0:nii.hdr.dime.dim(2)-1;
    ys = 0:nii.hdr.dime.dim(3)-1;
    zs = 0:nii.hdr.dime.dim(4)-1;
else
    [xs, ys, zs] = coords_from_nii(nii);
    if reorientate
        xs = -xs;
        ys = -ys;
    end
    %check for single slices
    if length(xs) == 1
        xs = xs + [-1 1]*nii.hdr.dime.pixdim(2)/6;
    end
    if length(ys) == 1
        ys = ys + [-1 1]*nii.hdr.dime.pixdim(3)/6;
    end
    if length(zs) == 1
        zs = zs + [-1 1]*nii.hdr.dime.pixdim(4)/6;
    end
end


%select slice and set coords for slice
switch lower(dim)
    case 'x'
        im = squeeze(V(slice,:,:))';
        xs_im = ys;
        ys_im = zs;
    case 'y'
        im = squeeze(V(:,slice,:))';
        xs_im = xs;
        ys_im = zs;
    case 'z'
        im = V(:,:,slice)';
        xs_im = xs;
        ys_im = ys;
    otherwise
        error('axis must be x, y, or z');
end

if switch_dims
    tmp = xs_im;
    xs_im = ys_im;
    ys_im = tmp;
    im = im';
end

%display slice and set colormap
h = imagesc(xs_im,ys_im,im,clim);
if ~show_axes
    ax=gca();
    ax.Visible=0;
end
colormap(cmap);
%ax=gca();
%ax.Colormap=cmap; %nonexistent 'Colormap' field in MATLAB 2017

% %convert slice to an index image into cmap, setting max and min values
% %according to clim
% im = grayslice(mat2gray(im,double(clim)),length(cmap));
% 
% %display slice using subimage - this converts an index image into a
% %truecolor image using the supplied colormap - this allows images to be
% %displayed in the same figure with different colourmaps - useful for colour
% %overlays
% h = subimage(xs_im,ys_im,im,cmap);

%set axis directions
if ~isequal(dim,'z')
    set(gca,'ydir','normal');
end

%scale axis if using vox coords
axis image
if use_vox_coords
    switch lower(dim)
        case 'x'
            daspect([1./nii.hdr.dime.pixdim(3:4) 1]);
        case 'y'
            daspect([1./nii.hdr.dime.pixdim([2 4]) 1]);
        case 'z'
            daspect([1./nii.hdr.dime.pixdim(2:3) 1]);
    end
end