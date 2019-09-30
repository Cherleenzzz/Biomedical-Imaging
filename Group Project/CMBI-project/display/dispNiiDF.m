function dispNiiDF(df_nii,dim,slice,spacing,type,linespec,df_is_disp,reorientate)
%dispNiiDF(df_nii,dim,slice,spacing,type,linespec,df_is_disp)
%function to display a slice from a deformation field using either a 
%deformed grid or arrows - to be used in conjunction with dispNiiSlice
%
%INPUTS:
%   df_nii - nifti structure containing deformation/displacement field and
%       header info
%   dim - the dimension that the slice is taken from ('x', 'y', or 'z')
%   slice - the slice number to display
%   spacing - the spacing (in mm) between grid/arrows [5]
%   type - 'grid' or 'arrows' ['grid']
%   linespec - used for plotting grid/arrows ['g']
%   df_is_disp - if true df_nii contains a displacement field, if false it
%       contains a deformation field [false]
%   reorientate: if true use -ve of x and y coords - this displays 
%               images same as in NiftyView (default = true)
%
%Note: deformation field = coords + displacement field

%check for optional inputs
if ~exist('spacing','var') || isempty(spacing)
    spacing = 5;
end
if ~exist('type','var') || isempty(type)
    type = 'grid';
end
if ~exist('linespec','var') || isempty(linespec)
    linespec = 'g';
end
if ~exist('df_is_disp','var') || isempty(df_is_disp)
    df_is_disp = false;
end
if ~exist('reorientate','var') || isempty(reorientate)
    reorientate = true;
end

%check if 1 or 3 spacings
if length(spacing) == 1
    spacing = [spacing spacing spacing];
elseif length(spacing) ~= 3
    error('must specify none, one, or three spacing values');
end

%get vox coords
[xs, ys, zs] = coords_from_nii(df_nii);
if reorientate
    xs = -xs;
    ys = -ys;
    df_nii.img(:,:,:,:,1:2) = -df_nii.img(:,:,:,:,1:2);
end

%select slice and set coords and def field for slice
switch lower(dim)
    case 'x'
        df_X = squeeze(df_nii.img(slice,:,:,1,2));
        df_Y = squeeze(df_nii.img(slice,:,:,1,3));
        xs_im = ys;
        ys_im = zs;
        spacing_x = spacing(2);
        spacing_y = spacing(3);
    case 'y'
        df_X = squeeze(df_nii.img(:,slice,:,1,1));
        df_Y = squeeze(df_nii.img(:,slice,:,1,3));
        xs_im = xs;
        ys_im = zs;
        spacing_x = spacing(1);
        spacing_y = spacing(3);
    case 'z'
        df_X = squeeze(df_nii.img(:,:,slice,1,1));
        df_Y = squeeze(df_nii.img(:,:,slice,1,2));
        xs_im = xs;
        ys_im = ys;
        spacing_x = spacing(1);
        spacing_y = spacing(2);
    otherwise
        error('dim must be x, y, or z');
end

%display def field depending on type
lims = axis;
[X,Y] = ndgrid(xs_im,ys_im);
df_X = double(df_X);
df_Y = double(df_Y);
if df_is_disp
    df_X = df_X + X;
    df_Y = df_Y + Y;
end
validPos=~isnan(df_X);
SI_X = scatteredInterpolant(df_X(validPos),df_Y(validPos),X(validPos));
SI_Y = scatteredInterpolant(df_X(validPos),df_Y(validPos),Y(validPos));
grid_xs_small = min(xs_im(:)):abs(xs_im(2)-xs_im(1)):max(xs_im(:));
grid_xs_large = min(xs_im(:)):spacing_x:max(xs_im(:));
grid_ys_small = min(ys_im(:)):abs(ys_im(2)-ys_im(1)):max(ys_im(:));
grid_ys_large = min(ys_im(:)):spacing_y:max(ys_im(:));
switch type
    case 'grid'
        hold on
        for x = 1:length(grid_xs_large)
            grid_disp_x = SI_X({grid_xs_large(x),grid_ys_small});
            grid_disp_y = SI_Y({grid_xs_large(x),grid_ys_small});
            plot(grid_disp_x,grid_disp_y,linespec);
        end
        for y = 1:length(grid_ys_large)
            grid_disp_x = SI_X({grid_xs_small,grid_ys_large(y)});
            grid_disp_y = SI_Y({grid_xs_small,grid_ys_large(y)});
            plot(grid_disp_x,grid_disp_y,linespec);
        end
        hold off
    case 'arrows'
        [arrow_X, arrow_Y] = ndgrid(grid_xs_large,grid_ys_large);
        arrow_disp_x = SI_X(arrow_X,arrow_Y) - arrow_X;
        arrow_disp_y = SI_Y(arrow_X,arrow_Y) - arrow_Y;
        hold on
        quiver(arrow_X,arrow_Y,arrow_disp_x,arrow_disp_y,0,linespec);
        hold off
    otherwise
        error('type must be grid or arrows');
end

%set axis directions
if ~isequal(dim,'z')
    set(gca,'ydir','normal');
end

if ~isequal(lims,[0 1 0 1])
    axis(lims);
end