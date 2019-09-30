% addpath
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

% load surrogate signal
load('surrogate.mat')

% use first 100 images
% load registration results

SI_def1 = zeros(100,1);
AP_def1 = zeros(100,1);
SI_def2 = zeros(100,1);
AP_def2 = zeros(100,1);
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
    CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
    
    SI_def1(i) = CPG1_nii.img(44,38,1,1,2);
    AP_def1(i) = CPG1_nii.img(44,38,1,1,1);
    SI_def2(i) = CPG2_nii.img(44,38,1,1,2);
    AP_def2(i) = CPG2_nii.img(44,38,1,1,1);
end

S_test = sub_pix(1:100);

%% linear
T = 1000;
para = zeros(2,T);
index = zeros(1,100);

for t = 1:T
    for k = 1:100
        index(k) = randi(100);    
    end
    A = SI_def1(index);
    S = [S_test(index) ones(100,1)];
    para(:,t) = S\A;
end

figure
hold on
histogram(para(1,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para(1,:))-2*std(para(1,:));
range_sigma(2) = mean(para(1,:))+2*std(para(1,:));
in_order = sort(para(1,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model1 c0')
hold off
sprintf('c0 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c0 95 range: %.3f to %.3f', range_95(1),range_95(2))

figure
hold on
histogram(para(2,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para(2,:))-2*std(para(2,:));
range_sigma(2) = mean(para(2,:))+2*std(para(2,:));
in_order = sort(para(2,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model1 c1')
hold off
sprintf('c1 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c1 95 range: %.3f to %.3f', range_95(1),range_95(2))




%% 2nd order 
para_2nd = zeros(3,T);  

for t = 1:T
    for k = 1:100
        index(k) = randi(100);    
    end
    A = SI_def1(index);
    S = [S_test(index).^2 S_test(index) ones(100,1)];
    para_2nd(:,t) = S\A;
end

figure       
hold on
histogram(para_2nd(1,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_2nd(1,:))-2*std(para_2nd(1,:));
range_sigma(2) = mean(para_2nd(1,:))+2*std(para_2nd(1,:));
in_order = sort(para_2nd(1,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model2 c0')
hold off
sprintf('c0 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c0 95 range: %.3f to %.3f', range_95(1),range_95(2))


figure       
hold on
histogram(para_2nd(2,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_2nd(2,:))-2*std(para_2nd(2,:));
range_sigma(2) = mean(para_2nd(2,:))+2*std(para_2nd(2,:));
in_order = sort(para_2nd(2,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model2 c1')
hold off
sprintf('c1 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c1 95 range: %.3f to %.3f', range_95(1),range_95(2))


figure       
hold on
histogram(para_2nd(3,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_2nd(3,:))-2*std(para_2nd(3,:));
range_sigma(2) = mean(para_2nd(3,:))+2*std(para_2nd(3,:));
in_order = sort(para_2nd(3,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model2 c2')
hold off
sprintf('c2 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c2 95 range: %.3f to %.3f', range_95(1),range_95(2))

%% 3rd order 
para_3rd = zeros(4,T);  

for t = 1:T
    for k = 1:100
        index(k) = randi(100);    
    end
    A = SI_def1(index);
    S = [S_test(index).^3 S_test(index).^2 S_test(index) ones(100,1)];
    para_3rd(:,t) = S\A;
end

figure       
hold on
histogram(para_3rd(1,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_3rd(1,:))-2*std(para_3rd(1,:));
range_sigma(2) = mean(para_3rd(1,:))+2*std(para_3rd(1,:));
in_order = sort(para_3rd(1,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model3 c0')
hold off
sprintf('c0 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c0 95 range: %.3f to %.3f', range_95(1),range_95(2))



figure       
hold on
histogram(para_3rd(2,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_3rd(2,:))-2*std(para_3rd(2,:));
range_sigma(2) = mean(para_3rd(2,:))+2*std(para_3rd(2,:));
in_order = sort(para_3rd(2,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model3 c1')
hold off
sprintf('c1 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c1 95 range: %.3f to %.3f', range_95(1),range_95(2))



figure       
hold on
histogram(para_3rd(3,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_3rd(3,:))-2*std(para_3rd(3,:));
range_sigma(2) = mean(para_3rd(3,:))+2*std(para_3rd(3,:));
in_order = sort(para_3rd(3,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model3 c2')
hold off
sprintf('c2 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c2 95 range: %.3f to %.3f', range_95(1),range_95(2))



figure       
hold on
histogram(para_3rd(4,:),20);
range_sigma = zeros(2,1);
range_95 = zeros(2,1);
range_sigma(1) = mean(para_3rd(4,:))-2*std(para_3rd(4,:));
range_sigma(2) = mean(para_3rd(4,:))+2*std(para_3rd(4,:));
in_order = sort(para_3rd(4,:));
range_95(1) = in_order(2.5*0.01*T);
range_95(2) = in_order(T-2.5*0.01*T);
plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
plot([range_95(1),range_95(2)],[165,165],'+-')
legend('P(x|A)','2考 range','95% range');
title('Model3 c3')
hold off
sprintf('c3 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
sprintf('c3 95 range: %.3f to %.3f', range_95(1),range_95(2))



% %% wild 
% T = 1000;
% para_wild = zeros(2,T);
% 
% S = [S_test ones(100,1)];
% c = S\SI_def1;
% E = SI_def1-S*c;
% 
% for t = 1:T   
%     A = S*c+E.*randn(100,1);
%     para_wild(:,t) = S\A;
% end
% 
% figure
% hold on
% histogram(para_wild(1,:),20);
% range_sigma = zeros(2,1);
% range_95 = zeros(2,1);
% range_sigma(1) = mean(para_wild(1,:))-2*std(para_wild(1,:));
% range_sigma(2) = mean(para_wild(1,:))+2*std(para_wild(1,:));
% in_order = sort(para_wild(1,:));
% range_95(1) = in_order(2.5*0.01*T);
% range_95(2) = in_order(T-2.5*0.01*T);
% plot([range_sigma(1),range_sigma(2)],[160,160],'+-')
% plot([range_95(1),range_95(2)],[165,165],'+-')
% legend('P(x|A)','2考 range','95% range');
% title('c0')
% hold off
% sprintf('c0 2 sigma range: %.3f to %.3f', range_sigma(1),range_sigma(2))
% sprintf('c0 95 range: %.3f to %.3f', range_95(1),range_95(2))


