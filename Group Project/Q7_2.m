% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

% load surrogate signal
load('surrogate.mat')

grad = gradient(sub_pix(2:101));


gradient_in = grad(grad<0);
gradient_ex = grad(grad>0);
inhalations = sub_pix(grad<0); 
exhalations = sub_pix(grad>0);

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

SI_def_in = SI_def(grad<0);
SI_def_ex = SI_def(grad>0);

n_in = length(inhalations);
n_ex = length(exhalations);

%% inhalations
% fit 1D linear model
S_lin_in = [inhalations ones(n_in, 1)];
c_lin_in = S_lin_in\SI_def_in;
y_lin_in = S_lin_in*c_lin_in;

SSD_lin_in = sum((y_lin_in-SI_def_in).^2);


% fit 2nd polynomial model
S_2nd_in = [inhalations.^2 inhalations ones(n_in,1)];
c_2nd_in = S_2nd_in\SI_def_in;
y_2nd_in = S_2nd_in*c_2nd_in;

SSD_2nd_in = sum((y_2nd_in-SI_def_in).^2);


% fit 3rd polynomial model
S_3rd_in = [inhalations.^3 inhalations.^2 inhalations ones(n_in,1)];
c_3rd_in = S_3rd_in\SI_def_in;
y_3rd_in = S_3rd_in*c_3rd_in;

SSD_3rd_in = sum((y_3rd_in-SI_def_in).^2);

%% exhalations
% fit 1D linear model
S_lin_ex = [exhalations ones(n_ex, 1)];
c_lin_ex = S_lin_ex\SI_def_ex;
y_lin_ex = S_lin_ex*c_lin_ex;

SSD_lin_ex = sum((y_lin_ex-SI_def_ex).^2);

% fit 2nd polynomial model
S_2nd_ex = [exhalations.^2 exhalations ones(n_ex,1)];
c_2nd_ex = S_2nd_ex\SI_def_ex;
y_2nd_ex = S_2nd_ex*c_2nd_ex;

SSD_2nd_ex = sum((y_2nd_ex-SI_def_ex).^2);

% fit 3rd polynomial model
S_3rd_ex = [exhalations.^3 exhalations.^2 exhalations ones(n_ex,1)];
c_3rd_ex = S_3rd_ex\SI_def_ex;
y_3rd_ex = S_3rd_ex*c_3rd_ex;

SSD_3rd_ex = sum((y_3rd_ex-SI_def_ex).^2);

SSD_lin=SSD_lin_in+SSD_lin_ex;
SSD_2nd=SSD_2nd_in+SSD_2nd_ex;
SSD_3rd=SSD_3rd_in+SSD_3rd_ex;

save('models8_10', 'c_lin_in','c_2nd_in','c_3rd_in','c_lin_ex','c_2nd_ex','c_3rd_ex');

figure();
plot(inhalations, y_lin_in,'b')
title('Model 8')
xlabel('surrogate signal');
ylabel('control-point value')
hold on
scatter(inhalations, SI_def_in, 40, 'b', 'Marker', 'o')
plot(exhalations, y_lin_ex, 'm')
scatter(exhalations, SI_def_ex, 45, 'm', 'Marker', 'x')
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'-ob');
h(2) = plot(NaN,NaN,'-xm');
legend(h, 'inhalation','exhalation');

