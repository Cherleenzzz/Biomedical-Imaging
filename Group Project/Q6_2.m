% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

N = 1500;
sub_pix = zeros(N,1);

for i=1:N
    if i<10
        f_name = sprintf('000%i.nii', i);
    elseif (10<=i) && (i<100) 
        f_name = sprintf('00%i.nii', i);
    elseif  (100<=i) && (i<1000)
        f_name = sprintf('0%i.nii', i);
    else 
        f_name = sprintf('%i.nii', i);
    end
nii = load_untouch_nii(f_name);

% extract values for the 1st z-slice: code taken from dispNiiSlice
V = nii.img;
intensity_y150 = V(:,150,1);

% find the first pixel with an intensity greater than 20
I = find(intensity_y150>=20);

% interpolate to subpixel accuracy
pix = I(1);
sub_pix(i) = (pix-1) + ((20 - intensity_y150(pix-1))/(intensity_y150(pix) - intensity_y150(pix-1)));
end

figure,
x = linspace(1,1500,1500);
plot(x,sub_pix)
xlim([0,300])
ylim([52.3,53.6])
xlabel('image number')
ylabel('skin position [index]')

source_nii = load_untouch_nii('0007.nii');
figure,
dispNiiSlice(source_nii,'z',1,[],[],[],[],true,[])
xlabel('AP')
ylabel('SI')
hold on
plot(I(1),150,'g*')

save('surrogate', 'sub_pix')


