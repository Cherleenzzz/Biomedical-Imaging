% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')
% load landmark file
load('landmark_pos_phys.mat')
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
[SI_reg1_def,AP_reg1_def,SI_reg2_def,AP_reg2_def] = deal(zeros(100,N_CP));

% all cp
for i=1:100
    SI_reg1_def(i,:)=reshape(CPG1_nii{i}.img(:,:,1,1,2),1,[]);
    AP_reg1_def(i,:)=reshape(CPG1_nii{i}.img(:,:,1,1,1),1,[]);
    SI_reg2_def(i,:)=reshape(CPG2_nii{i}.img(:,:,1,1,2),1,[]);
    AP_reg2_def(i,:)=reshape(CPG2_nii{i}.img(:,:,1,1,1),1,[]);
end
%% train
% linear model
S = [sub_pix(1:100) ones(100,1)];
c_lin_si1 = S\SI_reg1_def;
c_lin_ap1 = S\AP_reg1_def;

c_lin_si2 = S\SI_reg2_def;
c_lin_ap2 = S\AP_reg2_def;

% 2nd order polynomial
S = [sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_2nd_si1 = S\SI_reg1_def;
c_2nd_ap1 = S\AP_reg1_def;

c_2nd_si2 = S\SI_reg2_def;
c_2nd_ap2 = S\AP_reg2_def;

% 3rd order polynomial
S = [sub_pix(1:100).^3 sub_pix(1:100).^2 sub_pix(1:100) ones(100,1)];
c_3rd_si1 = S\SI_reg1_def;
c_3rd_ap1 = S\AP_reg1_def;

c_3rd_si2 = S\SI_reg2_def;
c_3rd_ap2 = S\AP_reg2_def;

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
% the control point grid as a nifti structure
CP1_lin = cell(1400,1);
CP1_2nd = cell(1400,1);
CP1_3rd = cell(1400,1);
CP2_lin = cell(1400,1);
CP2_2nd = cell(1400,1);
CP2_3rd = cell(1400,1);
for i=1:1400
    CP1_lin{i}=CPG1_nii{i+100};
    CP1_2nd{i}=CPG1_nii{i+100};
    CP1_3rd{i}=CPG1_nii{i+100};
    
    CP2_lin{i}=CPG2_nii{i+100};
    CP2_2nd{i}=CPG2_nii{i+100};
    CP2_3rd{i}=CPG2_nii{i+100};
end

for i = 1:1400
    % linear
    fit1_lin(:,:,1,1,1) = reshape(fit_lin_ap1(i,:), 67, 67, 1, 1, 1);
    fit1_lin(:,:,1,1,2) = reshape(fit_lin_si1(i,:), 67, 67, 1, 1, 1);
    CP1_lin{i}.img = fit1_lin;
    fit2_lin(:,:,1,1,1) = reshape(fit_lin_ap2(i,:), 67, 67, 1, 1, 1);
    fit2_lin(:,:,1,1,2) = reshape(fit_lin_si2(i,:), 67, 67, 1, 1, 1);
    CP2_lin{i}.img = fit2_lin;
    
    % 2nd
    fit1_2nd(:,:,1,1,1) = reshape(fit_2nd_ap1(i,:), 67, 67, 1, 1, 1);
    fit1_2nd(:,:,1,1,2) = reshape(fit_2nd_si1(i,:), 67, 67, 1, 1, 1);
    CP1_2nd{i}.img = fit1_2nd;
    fit2_2nd(:,:,1,1,1) = reshape(fit_2nd_ap2(i,:), 67, 67, 1, 1, 1);
    fit2_2nd(:,:,1,1,2) = reshape(fit_2nd_si2(i,:), 67, 67, 1, 1, 1);
    CP2_2nd{i}.img = fit2_2nd;
    
    % 3rd
    fit1_3rd(:,:,1,1,1) = reshape(fit_3rd_ap1(i,:), 67, 67, 1, 1, 1);
    fit1_3rd(:,:,1,1,2) = reshape(fit_3rd_si1(i,:), 67, 67, 1, 1, 1);
    CP1_3rd{i}.img = fit1_3rd;
    fit2_3rd(:,:,1,1,1) = reshape(fit_3rd_ap2(i,:), 67, 67, 1, 1, 1);
    fit2_3rd(:,:,1,1,2) = reshape(fit_3rd_si2(i,:), 67, 67, 1, 1, 1);
    CP2_3rd{i}.img = fit2_3rd;
end
% landmark estimation for all test images
est_landmark_lin = zeros(4,2,1400);
est_landmark_2nd = zeros(4,2,1400);
est_landmark_3rd = zeros(4,2,1400);
dist_nii = load_untouch_nii('0007_sdt.nii');
for i = 1:1400
    est_landmark_lin(:,:,i) = transPointsWithCPGsSliding(CP1_lin{i},CP2_lin{i},dist_nii,landmark_pos_phys(:,:,i+100),Image_nii{i+100}, false, true);
    est_landmark_2nd(:,:,i) = transPointsWithCPGsSliding(CP1_2nd{i},CP2_2nd{i},dist_nii,landmark_pos_phys(:,:,i+100),Image_nii{i+100}, false, true);
    est_landmark_3rd(:,:,i) = transPointsWithCPGsSliding(CP1_3rd{i},CP2_3rd{i},dist_nii,landmark_pos_phys(:,:,i+100),Image_nii{i+100}, false, true);
end

%% error for landmark1
L1_lin_si_e=sum(abs(landmark_pos_phys(1,2,101:1500)-est_landmark_lin(1,2,:)))/1400.0
L1_lin_ap_e=sum(abs(landmark_pos_phys(1,1,101:1500)-est_landmark_lin(1,1,:)))/1400.0
L1_quad_si_e=sum(abs(landmark_pos_phys(1,2,101:1500)-est_landmark_2nd(1,2,:)))/1400.0
L1_quad_ap_e=sum(abs(landmark_pos_phys(1,1,101:1500)-est_landmark_2nd(1,1,:)))/1400.0
L1_cubic_si_e=sum(abs(landmark_pos_phys(1,2,101:1500)-est_landmark_3rd(1,2,:)))/1400.0
L1_cubic_ap_e=sum(abs(landmark_pos_phys(1,1,101:1500)-est_landmark_3rd(1,1,:)))/1400.0

a=landmark_pos_phys(1,:,101:1500)-est_landmark_lin(1,:,:);
b=reshape(a,2,1400);
L1_lin_error=sqrt(sum(b.^2));
L1_norm2_lin=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(1,:,101:1500)-est_landmark_2nd(1,:,:);
b=reshape(a,2,1400);
L1_2nd_error=sqrt(sum(b.^2));
L1_norm2_2nd=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(1,:,101:1500)-est_landmark_3rd(1,:,:);
b=reshape(a,2,1400);
L1_3rd_error=sqrt(sum(b.^2));
L1_norm2_3rd=sum(sqrt(sum(b.^2)))/1400

%% error for landmark2
L2_lin_si_e=sum(abs(landmark_pos_phys(2,2,101:1500)-est_landmark_lin(2,2,:)))/1400.0
L2_lin_ap_e=sum(abs(landmark_pos_phys(2,1,101:1500)-est_landmark_lin(2,1,:)))/1400.0
L2_quad_si_e=sum(abs(landmark_pos_phys(2,2,101:1500)-est_landmark_2nd(2,2,:)))/1400.0
L2_quad_ap_e=sum(abs(landmark_pos_phys(2,1,101:1500)-est_landmark_2nd(2,1,:)))/1400.0
L2_cubic_si_e=sum(abs(landmark_pos_phys(2,2,101:1500)-est_landmark_3rd(2,2,:)))/1400.0
L2_cubic_ap_e=sum(abs(landmark_pos_phys(2,1,101:1500)-est_landmark_3rd(2,1,:)))/1400.0

a=landmark_pos_phys(2,:,101:1500)-est_landmark_lin(2,:,:);
b=reshape(a,2,1400);
L2_lin_error=sqrt(sum(b.^2));
L2_norm2_lin=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(2,:,101:1500)-est_landmark_2nd(2,:,:);
b=reshape(a,2,1400);
L2_2nd_error=sqrt(sum(b.^2));
L2_norm2_2nd=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(2,:,101:1500)-est_landmark_3rd(2,:,:);
b=reshape(a,2,1400);
L2_3rd_error=sqrt(sum(b.^2));
L2_norm2_3rd=sum(sqrt(sum(b.^2)))/1400

%% error for landmark3
L3_lin_si_e=sum(abs(landmark_pos_phys(3,2,101:1500)-est_landmark_lin(3,2,:)))/1400.0
L3_lin_ap_e=sum(abs(landmark_pos_phys(3,1,101:1500)-est_landmark_lin(3,1,:)))/1400.0
L3_quad_si_e=sum(abs(landmark_pos_phys(3,2,101:1500)-est_landmark_2nd(3,2,:)))/1400.0
L3_quad_ap_e=sum(abs(landmark_pos_phys(3,1,101:1500)-est_landmark_2nd(3,1,:)))/1400.0
L3_cubic_si_e=sum(abs(landmark_pos_phys(3,2,101:1500)-est_landmark_3rd(3,2,:)))/1400.0
L3_cubic_ap_e=sum(abs(landmark_pos_phys(3,1,101:1500)-est_landmark_3rd(3,1,:)))/1400.0

a=landmark_pos_phys(3,:,101:1500)-est_landmark_lin(3,:,:);
b=reshape(a,2,1400);
L3_lin_error=sqrt(sum(b.^2));
L3_norm2_lin=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(3,:,101:1500)-est_landmark_2nd(3,:,:);
b=reshape(a,2,1400);
L3_2nd_error=sqrt(sum(b.^2));
L3_norm2_2nd=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(3,:,101:1500)-est_landmark_3rd(3,:,:);
b=reshape(a,2,1400);
L3_3rd_error=sqrt(sum(b.^2));
L3_norm2_3rd=sum(sqrt(sum(b.^2)))/1400

%% error for landmark4
L4_lin_si_e=sum(abs(landmark_pos_phys(4,2,101:1500)-est_landmark_lin(4,2,:)))/1400.0
L4_lin_ap_e=sum(abs(landmark_pos_phys(4,1,101:1500)-est_landmark_lin(4,1,:)))/1400.0
L4_quad_si_e=sum(abs(landmark_pos_phys(4,2,101:1500)-est_landmark_2nd(4,2,:)))/1400.0
L4_quad_ap_e=sum(abs(landmark_pos_phys(4,1,101:1500)-est_landmark_2nd(4,1,:)))/1400.0
L4_cubic_si_e=sum(abs(landmark_pos_phys(4,2,101:1500)-est_landmark_3rd(4,2,:)))/1400.0
L4_cubic_ap_e=sum(abs(landmark_pos_phys(4,1,101:1500)-est_landmark_3rd(4,1,:)))/1400.0

a=landmark_pos_phys(4,:,101:1500)-est_landmark_lin(4,:,:);
b=reshape(a,2,1400);
L4_lin_error=sqrt(sum(b.^2));
L4_norm2_lin=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(4,:,101:1500)-est_landmark_2nd(4,:,:);
b=reshape(a,2,1400);
L4_2nd_error=sqrt(sum(b.^2));
L4_norm2_2nd=sum(sqrt(sum(b.^2)))/1400

a=landmark_pos_phys(4,:,101:1500)-est_landmark_3rd(4,:,:);
b=reshape(a,2,1400);
L4_3rd_error=sqrt(sum(b.^2));
L4_norm2_3rd=sum(sqrt(sum(b.^2)))/1400

% figure,
% plot(L1_lin_error,'r')
% title('Landmark1 error')
% hold on
% plot(L1_2nd_error,'b')
% hold on
% plot(L1_3rd_error,'g')
% xlim([600,800])
% legend('liner','2nd','3rd')
% hold off
% 
% figure,
% plot(L2_lin_error,'r')
% title('Landmark2 error')
% hold on
% plot(L2_2nd_error,'b')
% hold on
% plot(L2_3rd_error,'g')
% xlim([600,800])
% legend('liner','2nd','3rd')
% hold off
% 
% figure,
% plot(L3_lin_error,'r')
% title('Landmark3 error')
% hold on
% plot(L3_2nd_error,'b')
% hold on
% plot(L3_3rd_error,'g')
% xlim([600,800])
% legend('liner','2nd','3rd')
% hold off
% 
% figure,
% plot(L4_lin_error,'r')
% title('Landmark4 error')
% hold on
% plot(L4_2nd_error,'b')
% hold on
% plot(L4_3rd_error,'g')
% xlim([600,800])
% legend('liner','2nd','3rd')
% hold off

F(50) = struct('cdata',[],'colormap',[]);
figure,
set (gcf,'Position',[200 200 1000 500]);
for i = 701 : 750
    % load image
    current_image = Image_nii{100+i};
    
    subplot(1,3,1)
    dispNiiSlice(current_image,'z',1,[],[],[],[],false,[])
    title('Estimated liner model landmarks')
    hold on
    plot(landmark_pos_phys(:,1,i), landmark_pos_phys(:,2,i),'r*', 'markersize', 1)
    hold on
    plot(est_landmark_lin(:,1,i), est_landmark_lin(:,2,i),'b*', 'markersize', 1)
    legend('truth landmarks','estimate liner landmarks')
    hold off
    
    subplot(1,3,2)
    dispNiiSlice(current_image,'z',1,[],[],[],[],false,[])
    title('Estimated 1D 2nd model landmarks')
    hold on
    plot(landmark_pos_phys(:,1,i), landmark_pos_phys(:,2,i),'r*', 'markersize', 1)
    hold on
    plot(est_landmark_2nd(:,1,i), est_landmark_2nd(:,2,i),'b*', 'markersize', 1)
    legend('truth landmarks','estimate 1D 2nd landmarks')
    hold off
    
    subplot(1,3,3);
    dispNiiSlice(current_image,'z',1,[],[],[],[],false,[])
    title('Estimated 1D 3rd model landmarks')
    hold on
    plot(landmark_pos_phys(:,1,i), landmark_pos_phys(:,2,i),'r*', 'markersize', 1)
    hold on
    plot(est_landmark_3rd(:,1,i), est_landmark_3rd(:,2,i),'b*', 'markersize', 1)
    legend('truth landmarks','estimate 1D 3rd landmarks')
    hold off
   
    % save frames
    F(i) = getframe(gcf);
    
    pause(0.01)
    
end
SaveAsVideo(F,'landmark.mp4');