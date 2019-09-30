% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

N = 1500;
sub_pix1 = zeros(N,1);
sub_pix2 = zeros(N,1);

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
intensity_y75 = V(:,75,1);

% find the first pixel with an intensity greater than 20
I1 = find(intensity_y150>=20);
I2 = find(intensity_y75>=20);

% interpolate to subpixel accuracy
pix1 = I1(1);
sub_pix1(i) = (pix1-1) + ((20 - intensity_y150(pix1-1))/(intensity_y150(pix1) - intensity_y150(pix1-1)));

pix2 = I2(1);
sub_pix2(i) = (pix2-1) + ((20 - intensity_y75(pix2-1))/(intensity_y75(pix2) - intensity_y75(pix2-1)));
end

x = linspace(1,1500,1500);
plot(x,sub_pix1,'r-')
hold on
plot(x,sub_pix2, 'g-')
xlim([0,300])
xlabel('image number')
ylabel('skin position [index]')
legend('y=150', 'y=75')

save('surrogate1', 'sub_pix1')
save('surrogate2', 'sub_pix2')

% %% fit models
% % use first 100 images
% % load registration results
% SI_def = zeros(100,1);
% for i=1:100
%     if i<10
%         f_name = sprintf('000%i', i);
%     elseif (10<=i) && (i<100)
%         f_name = sprintf('00%i', i);
%     elseif  (100<=i) && (i<1000)
%         f_name = sprintf('0%i', i);
%     else
%         f_name = sprintf('%i', i);
%     end
%     
%     CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
%     
%     SI_def(i) = CPG1_nii.img(44,38,1,1,2);
% end
% 
% % fit 1D linear model
% S = [sub_pix2(1:100) ones(100,1)];
% c = S\SI_def;
% y = S*c;
% 
% % fit 2nd polynomial model
% S_2nd = [sub_pix2(1:100).^2 sub_pix2(1:100) ones(100,1)];
% c_2nd = S_2nd\SI_def;
% y_2nd = S_2nd*c_2nd;
% 
% % fit 3rd polynomial model
% S_3rd = [sub_pix2(1:100).^3 sub_pix2(1:100).^2 sub_pix2(1:100) ones(100,1)];
% c_3rd = S_3rd\SI_def;
% y_3rd = S_3rd*c_3rd;
% 
% % plot CP deformation against surrogate signal
% figure
% plot(sub_pix2(1:100), SI_def, '-x')
% xlabel('surrogate signal')
% ylabel('control-point value')
% hold on
% 
% % plot 1D linear model fit
% plot(sub_pix2(1:100),y,'r')
% hold on
% % plot 2nd polynomial model fit
% plot(sub_pix2(1:100),y_2nd,'k')
% hold on
% % plot 3rd polynomial model fit
% plot(sub_pix2(1:100),y_3rd,'g')
% legend('CP deformation','1D liner model','2nd polynomial model','3rd polynomial model');
% 
% % residual fitting error
% SSD_lin = sum((y-SI_def).^2);
% SSD_2nd = sum((y_2nd-SI_def).^2);
% SSD_3rd = sum((y_3rd-SI_def).^2);
% 
% 
% %% fit model for every CP deformation
% 
% % number of CP
% N_CP = 67^2;
% 
% % initialise matricies
% [SI_reg1_def,AP_reg1_def,SI_reg2_def,AP_reg2_def] = deal(zeros(100,N_CP));
% for i=1:100
%     if i<10
%         f_name = sprintf('000%i', i);
%     elseif (10<=i) && (i<100)
%         f_name = sprintf('00%i', i);
%     elseif (100<=i) && (i<1000)
%         f_name = sprintf('0%i', i);
%     else
%         f_name = sprintf('%i', i);
%     end
%     CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
%     CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
%     SI_reg1_def(i,:)=reshape(CPG1_nii.img(:,:,1,1,2),1,[]);
%     AP_reg1_def(i,:)=reshape(CPG1_nii.img(:,:,1,1,1),1,[]);
%     SI_reg2_def(i,:)=reshape(CPG2_nii.img(:,:,1,1,2),1,[]);
%     AP_reg2_def(i,:)=reshape(CPG2_nii.img(:,:,1,1,1),1,[]);
% end
% 
% X = [SI_reg1_def AP_reg1_def SI_reg2_def AP_reg2_def];
% 
% % linear model
% S = [sub_pix2(1:100) ones(100,1)];
% c_lin = S\X;
% y_lin = S*c_lin;
% 
% SSD_lin_all = sum((y_lin-X).^2);
% 
% % 2nd order polynomial
% S = [sub_pix2(1:100).^2 sub_pix2(1:100) ones(100,1)];
% c_2nd = S\X;
% y_2nd = S*c_2nd;
% 
% SSD_2nd_all = sum((y_2nd-X).^2);
% 
% % 3rd order polynomial
% S = [sub_pix2(1:100).^3 sub_pix2(1:100).^2 sub_pix2(1:100) ones(100,1)];
% c_3rd = S\X;
% y_3rd = S*c_3rd;
% 
% SSD_3rd_all = sum((y_3rd-X).^2);
% 
% save('models', 'c_lin','c_2nd','c_3rd');

%% fit models
% use first 100 images
% load registration results
grad = gradient(sub_pix1(2:101));
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

% fit linear polynomial model
S_lin4 = [sub_pix2(1:100) sub_pix1(1:100) ones(100,1)];
c_lin4 = S_lin4\SI_def;
y_lin4 = S_lin4*c_lin4;

S_lin6 = [grad sub_pix1(1:100) ones(100,1)];
c_lin6 = S_lin6\SI_def;
y_lin6 = S_lin6*c_lin6;

S_2nd5 = [sub_pix2(1:100).^2 sub_pix1(1:100) ones(100,1)];
c_2nd5 = S_2nd5\SI_def;
y_2nd5 = S_2nd5*c_2nd5;

S_2nd7 = [grad.^2 sub_pix1(1:100) ones(100,1)];
c_2nd7 = S_2nd7\SI_def;
y_2nd7 = S_2nd7*c_2nd7;

SSD_4 = sum((y_lin4-SI_def).^2);
SSD_5 = sum((y_2nd5-SI_def).^2);
SSD_6 = sum((y_lin6-SI_def).^2);
SSD_7 = sum((y_2nd7-SI_def).^2);



% plot CP deformation against surrogate signal
figure
plot3(sub_pix1(1:100),grad, SI_def, '-r*')
xlabel('surrogate signal 1')
ylabel('gradient')
zlabel('control-point value')
title ('Model 7')
hold on

xx=linspace(min(sub_pix1(1:100)),max(sub_pix1(1:100)));

yy=linspace(min(grad),max(grad));

[x,y]=meshgrid(xx,yy);

z=(-14.957824821903118)*x+3.046152855271652*(y.^2)+(8.064620590174280e+02)*ones(100,100);

surf(x,y,z);
grid on;

% % use first 100 images
% % load registration results
% grad = gradient(sub_pix1(2:101));
% SI_def = zeros(100,1);
% for i=1:100
%     if i<10
%         f_name = sprintf('000%i', i);
%     elseif (10<=i) && (i<100)
%         f_name = sprintf('00%i', i);
%     elseif  (100<=i) && (i<1000)
%         f_name = sprintf('0%i', i);
%     else
%         f_name = sprintf('%i', i);
%     end
%     
%     CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
%     
%     SI_def(i) = CPG1_nii.img(44,38,1,1,2);
% end
% 
% 
% % fit 2nd polynomial model
% S_2nd = [grad.^2 sub_pix1(1:100) ones(100,1)];
% c_2nd = S_2nd\SI_def;
% y_2nd = S_2nd*c_2nd;


% % residual fitting error
% SSD_2nd = sum((y_2nd-SI_def).^2);



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

% 2nd order polynomial
S_lin4 = [sub_pix2(1:100) sub_pix1(1:100) ones(100,1)];
c_lin4 = S_lin4\X;
y_lin4 = S_lin4*c_lin4;

S_lin6 = [grad sub_pix1(1:100) ones(100,1)];
c_lin6 = S_lin6\X;
y_lin6 = S_lin6*c_lin6;

S_2nd5 = [sub_pix2(1:100).^2 sub_pix1(1:100) ones(100,1)];
c_2nd5 = S_2nd5\X;
y_2nd5 = S_2nd5*c_2nd5;

S_2nd7 = [grad.^2 sub_pix1(1:100) ones(100,1)];
c_2nd7 = S_2nd7\X;
y_2nd7 = S_2nd7*c_2nd7;

save('models4_7', 'c_lin4','c_2nd5','c_lin6','c_2nd7');
SSD_4_all = sum(sum((y_lin4-X).^2));
SSD_5_all = sum(sum((y_2nd5-X).^2));
SSD_6_all = sum(sum((y_lin6-X).^2));
SSD_7_all = sum(sum((y_2nd7-X).^2));