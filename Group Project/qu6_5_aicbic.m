% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')
% load surrogate signal
load('surrogate.mat')

% load image_nii
Image_nii = cell(1500,1);
for i=1:1500
    if i<10
        f_name = sprintf('000%i', i);
    elseif (10<=i) && (i<100)
        f_name = sprintf('00%i', i);
    elseif  (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else
        f_name = sprintf('%i', i);
    end
    Image_nii{i} = load_untouch_nii(sprintf('%s.nii', f_name));
end

% load CP_nii
CPG1_nii = cell(1500,1);
CPG2_nii = cell(1500,1);
for i=1:1500
    if i<10
        f_name = sprintf('000%i', i);
    elseif (10<=i) && (i<100)
        f_name = sprintf('00%i', i);
    elseif (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else
        f_name = sprintf('%i', i);
    end
    CPG1_nii{i} = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    CPG2_nii{i} = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
end

% number of CP
N_CP = 67^2;

% initialise matricies
[SI_reg1_def_train,AP_reg1_def_train,SI_reg2_def_train,AP_reg2_def_train] = deal(zeros(100,N_CP));
[SI_reg1_def_test,AP_reg1_def_test,SI_reg2_def_test,AP_reg2_def_test] = deal(zeros(1400,N_CP));

% all cp
% train
for i=1:100
    SI_reg1_def_train(i,:)=reshape(CPG1_nii{i}.img(:,:,1,1,2),1,[]);
    AP_reg1_def_train(i,:)=reshape(CPG1_nii{i}.img(:,:,1,1,1),1,[]);
    SI_reg2_def_train(i,:)=reshape(CPG2_nii{i}.img(:,:,1,1,2),1,[]);
    AP_reg2_def_train(i,:)=reshape(CPG2_nii{i}.img(:,:,1,1,1),1,[]);
end
for i=1:1400
    SI_reg1_def_test(i,:)=reshape(CPG1_nii{i+100}.img(:,:,1,1,2),1,[]);
    AP_reg1_def_test(i,:)=reshape(CPG1_nii{i+100}.img(:,:,1,1,1),1,[]);
    SI_reg2_def_test(i,:)=reshape(CPG2_nii{i+100}.img(:,:,1,1,2),1,[]);
    AP_reg2_def_test(i,:)=reshape(CPG2_nii{i+100}.img(:,:,1,1,1),1,[]);
end
%% train
% linear model
S = [sub_pix(1:100) ones(100,1)];
c_lin_si1 = S\SI_reg1_def_train;
c_lin_ap1 = S\AP_reg1_def_train;

c_lin_si2 = S\SI_reg2_def_train;
c_lin_ap2 = S\AP_reg2_def_train;

% 2nd order polynomial
S = [sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_2nd_si1 = S\SI_reg1_def_train;
c_2nd_ap1 = S\AP_reg1_def_train;

c_2nd_si2 = S\SI_reg2_def_train;
c_2nd_ap2 = S\AP_reg2_def_train;

% 3rd order polynomial
S = [sub_pix(1:100).^3 sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_3rd_si1 = S\SI_reg1_def_train;
c_3rd_ap1 = S\AP_reg1_def_train;

c_3rd_si2 = S\SI_reg2_def_train;
c_3rd_ap2 = S\AP_reg2_def_train;

%% test
% linear model
S = [sub_pix(101:1500) ones(1400,1)];
fit_lin_si1 = S*c_lin_si1;
fit_lin_ap1 = S*c_lin_ap1;

fit_lin_si2 = S*c_lin_si2;
fit_lin_ap2 = S*c_lin_ap2;

% 2nd order polynomial
S = [sub_pix(101:1500).^2 sub_pix(101:1500) ones(1400,1)];
fit_2nd_si1 = S*c_2nd_si1;
fit_2nd_ap1 = S*c_2nd_ap1;

fit_2nd_si2 = S*c_2nd_si2;
fit_2nd_ap2 = S*c_2nd_ap2;

% 3rd order polynomial
S = [sub_pix(101:1500).^3 sub_pix(101:1500).^2 sub_pix(101:1500) ones(1400,1)];
fit_3rd_si1 = S*c_3rd_si1;
fit_3rd_ap1 = S*c_3rd_ap1;

fit_3rd_si2 = S*c_3rd_si2;
fit_3rd_ap2 = S*c_3rd_ap2;

% Compute RSS
% region 1
fit_r1_si = {fit_lin_si1; fit_2nd_si1; fit_3rd_si1};
rss_arr_r1_si = cell(3,1);

for i=1:3
    rss_arr_r1_si{i} = sum((SI_reg1_def_test - fit_r1_si{i}).^2, 1);
end

fit_r1_ap = {fit_lin_ap1; fit_2nd_ap1; fit_3rd_ap1};
rss_arr_r1_ap = cell(3,1);

for i=1:3
    rss_arr_r1_ap{i} = sum((AP_reg1_def_test - fit_r1_ap{i}).^2, 1);
end
% region 2
fit_r2_si = {fit_lin_si2; fit_2nd_si2; fit_3rd_si2};
rss_arr_r2_si = cell(3,1);

for i=1:3
    rss_arr_r2_si{i} = sum((SI_reg2_def_test - fit_r2_si{i}).^2, 1);
end

fit_r2_ap = {fit_lin_ap2; fit_2nd_ap2; fit_3rd_ap2};
rss_arr_r2_ap = cell(3,1);

for i=1:3
    rss_arr_r2_ap{i} = sum((AP_reg2_def_test - fit_r2_ap{i}).^2, 1);
end

% Compute AIC/BIC
% N: number of model parameters
% K: number of data points
N_para = {3; 4; 5};
K = 1400;

% region1
AIC_r1_si = zeros(67*67, 3);
BIC_r1_si = zeros(67*67, 3);

for i=1:3
    N = N_para{i};
    RSS = rss_arr_r1_si{i};
    AIC_r1_si(:,i) = 2*N + K.*log(RSS./K);
    BIC_r1_si(:,i) = N*log(K) + K.*log(RSS./K);
end
AIC_r1_si=sum(AIC_r1_si)/4489
BIC_r1_si=sum(BIC_r1_si)/4489

AIC_r1_ap = zeros(67*67, 3);
BIC_r1_ap = zeros(67*67, 3);

for i=1:3
    N = N_para{i};
    RSS = rss_arr_r1_ap{i};
    AIC_r1_ap(:,i) = 2*N + K.*log(RSS./K);
    BIC_r1_ap(:,i) = N*log(K) + K.*log(RSS./K);
end
AIC_r1_ap=sum(AIC_r1_ap)/4489
BIC_r1_ap=sum(BIC_r1_ap)/4489

% region2
AIC_r2_si = zeros(67*67, 3);
BIC_r2_si = zeros(67*67, 3);

for i=1:3
    N = N_para{i};
    RSS = rss_arr_r2_si{i};
    AIC_r2_si(:,i) = 2*N + K.*log(RSS./K);
    BIC_r2_si(:,i) = N*log(K) + K.*log(RSS./K);
end
AIC_r2_si=sum(AIC_r2_si)/4489
BIC_r2_si=sum(BIC_r2_si)/4489

AIC_r2_ap = zeros(67*67, 3);
BIC_r2_ap = zeros(67*67, 3);

for i=1:3
    N = N_para{i};
    RSS = rss_arr_r2_ap{i};
    AIC_r2_ap(:,i) = 2*N + K.*log(RSS./K);
    BIC_r2_ap(:,i) = N*log(K) + K.*log(RSS./K);
end
AIC_r2_ap=sum(AIC_r2_ap)/4489
BIC_r2_ap=sum(BIC_r2_ap)/4489