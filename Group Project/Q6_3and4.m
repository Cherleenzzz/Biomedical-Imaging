% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

% load surrogate signal
load('surrogate.mat')

% use first 100 images
% load registration results
SI_def = zeros(100,1);
for i=1:100
    if i<10
        f_name = sprintf('000%i', i);
    elseif (10<=i) && (i<100)
        f_name = sprintf('00%i', i);
    elseif  (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else
        f_name = sprintf('%i', i);
    end
    
    CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    
    SI_def(i) = CPG1_nii.img(44,38,1,1,2);
end

% fit 1D linear model
S = [sub_pix(1:100) ones(100,1)];
c = S\SI_def;
y = S*c;

% fit 2nd polynomial model
S_2nd = [sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_2nd = S_2nd\SI_def;
y_2nd = S_2nd*c_2nd;

% fit 3rd polynomial model
S_3rd = [sub_pix(1:100).^3 sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_3rd = S_3rd\SI_def;
y_3rd = S_3rd*c_3rd;

% plot CP deformation against surrogate signal
figure
plot(sub_pix(1:100), SI_def, '-x')
xlabel('surrogate signal')
ylabel('control-point value')
hold on

% plot 1D linear model fit
plot(sub_pix(1:100),y,'m')
title('Model 1, 2 and 3')
hold on
% plot 2nd polynomial model fit
plot(sub_pix(1:100),y_2nd,'k')
hold on
% plot 3rd polynomial model fit
plot(sub_pix(1:100),y_3rd,'g')
legend('CP deformation','1D liner model','2nd polynomial model','3rd polynomial model');

% residual fitting error
SSD_lin = sum((y-SI_def).^2);
SSD_2nd = sum((y_2nd-SI_def).^2);
SSD_3rd = sum((y_3rd-SI_def).^2);


%% fit model for every CP deformation

% number of CP
N_CP = 67^2;

% initialise matricies
[SI_reg1_def,AP_reg1_def,SI_reg2_def,AP_reg2_def] = deal(zeros(100,N_CP));
for i=1:100
    if i<10
        f_name = sprintf('000%i', i);
    elseif (10<=i) && (i<100)
        f_name = sprintf('00%i', i);
    elseif (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else
        f_name = sprintf('%i', i);
    end
    CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
    SI_reg1_def(i,:)=reshape(CPG1_nii.img(:,:,1,1,2),1,[]);
    AP_reg1_def(i,:)=reshape(CPG1_nii.img(:,:,1,1,1),1,[]);
    SI_reg2_def(i,:)=reshape(CPG2_nii.img(:,:,1,1,2),1,[]);
    AP_reg2_def(i,:)=reshape(CPG2_nii.img(:,:,1,1,1),1,[]);
end

X = [SI_reg1_def AP_reg1_def SI_reg2_def AP_reg2_def];

% linear model
S = [sub_pix(1:100) ones(100,1)];
c_lin = S\X;
y_lin = S*c_lin;

SSD_lin_all = sum(sum((y_lin-X).^2));

% 2nd order polynomial
S = [sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_2nd = S\X;
y_2nd = S*c_2nd;

SSD_2nd_all = sum(sum((y_2nd-X).^2));

% 3rd order polynomial
S = [sub_pix(1:100).^3 sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_3rd = S\X;
y_3rd = S*c_3rd;

SSD_3rd_all = sum(sum((y_3rd-X).^2));

save('models1_3', 'c_lin','c_2nd','c_3rd');





