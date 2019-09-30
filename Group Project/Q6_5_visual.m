% % add path
% addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
% addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
% addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
% addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
% addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')
% 
% 
% % N=1500-100=1400
% N=50;
% N_CP = 67^2;
% 
% 
% disp('Initialize - done')
% %% load datas
% % load surrogate signal as sub_pix
% load('surrogate.mat')
% % load models as c_1s, c_2nd, c_3rd
% load('models1_3.mat')
% 
% % load nii
% dist_nii = load_untouch_nii('0007_sdt.nii');
% source_nii = load_untouch_nii('0007.nii');
% 
% target_nii = [];
% CPG1_nii = [];
% CPG2_nii = [];
% 
% for i=101:100+N
%    if  (100<=i) && (i<1000)
%         f_name = sprintf('0%i', i);
%     else 
%         f_name = sprintf('%i', i);
%     end
%     target = load_untouch_nii(sprintf('%s.nii', f_name));
%     CPG1 = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
%     CPG2 = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
%     
%     target_nii=[target_nii;target];
%     CPG1_nii = [CPG1_nii,CPG1];
%     CPG2_nii = [CPG2_nii,CPG2];
% end
% 
% 
% disp('Load data - done')
% %% Calculate transformation
% S1 = [sub_pix(101:100+N) ones(N,1)];
% S2 = [sub_pix(101:100+N).^2 sub_pix(101:100+N) ones(N,1)];
% S3 = [sub_pix(101:100+N).^3 sub_pix(101:100+N).^2 sub_pix(101:100+N) ones(N,1)];
% x_1st = S1*c_lin;
% x_2nd = S2*c_2nd;
% x_3rd = S3*c_3rd;
% 
% disp('Calculate transformation - done')
% %% save into registrations
% [SI_reg1_1st,AP_reg1_1st,SI_reg2_1st,AP_reg2_1st] = deal(zeros(N,N_CP));
% [SI_reg1_2nd,AP_reg1_2nd,SI_reg2_2nd,AP_reg2_2nd] = deal(zeros(N,N_CP));
% [SI_reg1_3rd,AP_reg1_3rd,SI_reg2_3rd,AP_reg2_3rd] = deal(zeros(N,N_CP));
% 
% CPG1_nii_1st = CPG1_nii;
% CPG2_nii_1st = CPG2_nii;
% CPG1_nii_2nd = CPG1_nii;
% CPG2_nii_2nd = CPG2_nii;
% CPG1_nii_3rd = CPG1_nii;
% CPG2_nii_3rd = CPG2_nii;
% 
% for i=1:N
%     
%     % linear model
%     SI_reg1_1st(i,:) = x_1st(i,(1:N_CP));
%     AP_reg1_1st(i,:) = x_1st(i,(1+N_CP:2*N_CP));
%     SI_reg2_1st(i,:) = x_1st(i,(1+2*N_CP:3*N_CP));
%     AP_reg2_1st(i,:) = x_1st(i,(1+3*N_CP:4*N_CP));
%     CPG1_nii_1st(i).img(:,:,1,1,2) = reshape(SI_reg1_1st(i,:),[67,67]);
%     CPG1_nii_1st(i).img(:,:,1,1,1) = reshape(AP_reg1_1st(i,:),[67,67]);
%     CPG2_nii_1st(i).img(:,:,1,1,2) = reshape(SI_reg2_1st(i,:),[67,67]);
%     CPG2_nii_1st(i).img(:,:,1,1,1) = reshape(AP_reg2_1st(i,:),[67,67]);
%     
%     % second order model
%     SI_reg1_2nd(i,:) = x_2nd(i,(1:N_CP));
%     AP_reg1_2nd(i,:) = x_2nd(i,(1+N_CP:2*N_CP)); 
%     SI_reg2_2nd(i,:) = x_2nd(i,(1+2*N_CP:3*N_CP)); 
%     AP_reg2_2nd(i,:) = x_2nd(i,(1+3*N_CP:4*N_CP));
%     CPG1_nii_2nd(i).img(:,:,1,1,2) = reshape(SI_reg1_2nd(i,:),[67,67]);
%     CPG1_nii_2nd(i).img(:,:,1,1,1) = reshape(AP_reg1_2nd(i,:),[67,67]);
%     CPG2_nii_2nd(i).img(:,:,1,1,2) = reshape(SI_reg2_2nd(i,:),[67,67]);
%     CPG2_nii_2nd(i).img(:,:,1,1,1) = reshape(AP_reg2_2nd(i,:),[67,67]);
%     
%     % third order model
%     SI_reg1_3rd(i,:) = x_3rd(i,(1:N_CP));
%     AP_reg1_3rd(i,:) = x_3rd(i,(1+N_CP:2*N_CP));  
%     SI_reg2_3rd(i,:) = x_3rd(i,(1+2*N_CP:3*N_CP));  
%     AP_reg2_3rd(i,:) = x_3rd(i,(1+3*N_CP:4*N_CP));
%     
%     CPG1_nii_3rd(i).img(:,:,1,1,2) = reshape(SI_reg1_3rd(i,:),[67,67]);
%     CPG1_nii_3rd(i).img(:,:,1,1,1) = reshape(AP_reg1_3rd(i,:),[67,67]);
%     CPG2_nii_3rd(i).img(:,:,1,1,2) = reshape(SI_reg2_3rd(i,:),[67,67]);
%     CPG2_nii_3rd(i).img(:,:,1,1,1) = reshape(AP_reg2_3rd(i,:),[67,67]);
% end
% 
% disp('Save into registrations - done')
% %% present results & Residual motion
% F(N) = struct('cdata',[],'colormap',[]);
% Res_motion1 = zeros(N,160,160);
% Res_motion2 = zeros(N,160,160);
% Res_motion3 = zeros(N,160,160);
% 
% Def_error1=zeros(160,160,1,1,2,N);
% Def_error2=zeros(160,160,1,1,2,N);
% Def_error3=zeros(160,160,1,1,2,N);
% 
% for i=1:N
%     % registration     
%      [def_vol_nii, def_field_nii, dis_field_nii] = deformNiiWithCPGsSliding(CPG1_nii(i),CPG2_nii(i),...
%          dist_nii,source_nii,target_nii(i),false);
%      subplot(2,2,1)
%      dispNiiSlice(def_vol_nii,'z',1,[],[],[],[],true,[]);
%      title('Registration Image')
%     
%     % plot linear model
%     [def_vol_nii_1, def_field_nii_1, dis_field_nii_1] = deformNiiWithCPGsSliding(CPG1_nii_1st(i),CPG2_nii_1st(i),...
%         dist_nii,source_nii,target_nii(i),false);
%     subplot(2,2,2)
%     dispNiiSlice(def_vol_nii_1,'z',1,[],[],[],[],true,[]);
%     title('Model 1')
%     
%     % plot 2nd order model
%     [def_vol_nii_2, def_field_nii_2, dis_field_nii_2] = deformNiiWithCPGsSliding(CPG1_nii_2nd(i),CPG2_nii_2nd(i),...
%         dist_nii,source_nii,target_nii(i),false);
%     subplot(2,2,3)
%     dispNiiSlice(def_vol_nii_2,'z',1,[],[],[],[],true,[]);
%     title('Model 2')
%     
%     % plot 3rd order model
%     [def_vol_nii_3, def_field_nii_3, dis_field_nii_3] = deformNiiWithCPGsSliding(CPG1_nii_3rd(i),CPG2_nii_3rd(i),...
%         dist_nii,source_nii,target_nii(i),false);
%     subplot(2,2,4)
%     dispNiiSlice(def_vol_nii_3,'z',1,[],[],[],[],true,[]);
%     title('Model 3')
%     
%     % save frames
%     F(i) = getframe(gcf);
%     
%     % Residual motion
%     Res_motion1(i,:,:) = def_vol_nii_1.img-def_vol_nii.img;
%     Res_motion2(i,:,:) = def_vol_nii_2.img-def_vol_nii.img;
%     Res_motion3(i,:,:) = def_vol_nii_3.img-def_vol_nii.img;
% 
%     Def_error1(:,:,1,1,:,i)=def_field_nii_1.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
%     Def_error2(:,:,1,1,:,i)=def_field_nii_2.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
%     Def_error3(:,:,1,1,:,i)=def_field_nii_3.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
%     
%     
% end
% % SaveAsVideo(F,'myVideo.mp4');
% 
% 
% %% Residual motion Results
% mu_res_1 = zeros(N,1);
% mu_res_2 = zeros(N,1);
% mu_res_3 = zeros(N,1);
% for i=1:N
%      res_img_1 = reshape(Res_motion1(i,:,:),[160,160]);
%      res_img_2 = reshape(Res_motion2(i,:,:),[160,160]);
%      res_img_3 = reshape(Res_motion3(i,:,:),[160,160]);
%     subplot(1,3,1)
%     imshow(res_img_1');colormap('gray');title('Model 1');
%     subplot(1,3,2)
%     imshow(res_img_2');colormap('gray');title('Model 2');
%     subplot(1,3,3)
%     imshow(res_img_3');colormap('gray');title('Model 3');
%     mu_res_1(i) = nansum(nansum(abs(res_img_1)))/(160^2);
%     mu_res_2(i) = nansum(nansum(abs(res_img_2)))/(160^2);
%     mu_res_3(i) = nansum(nansum(abs(res_img_3)))/(160^2);
% end
% % figure
% % hold on
% % plot(mu_res_1);
% % plot(mu_res_2);
% % plot(mu_res_3);
% % legend('linear model','2nd order model','3rd order model')
% % title('Q6_5_1:residual motion')
% % saveas(gcf,'plots of residual motion.png')
% 
% save('Def_error1','Def_error1');
% save('Def_error2','Def_error2');
% save('Def_error3','Def_error3');
% 
% disp('Plots & Residual motion - done');
% 
% disp('####### Question 6.5 task 1 completed #######');

% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')


% N=1500-100=1400
N=50;
N_CP = 67^2;


disp('Initialize - done')
%% load datas
% load surrogate signal as sub_pix
load('surrogate1.mat')
load('surrogate2.mat')
% load models as c_1s, c_2nd, c_3rd
load('models4_7.mat')

% load nii
dist_nii = load_untouch_nii('0007_sdt.nii');
source_nii = load_untouch_nii('0007.nii');

target_nii = [];
CPG1_nii = [];
CPG2_nii = [];

%CPG2_nii = repmat(CPG2_nii,N,1);

for i=101:100+N
   if  (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else 
        f_name = sprintf('%i', i);
    end
    target = load_untouch_nii(sprintf('%s.nii', f_name));
    CPG1 = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    CPG2 = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
    
    target_nii=[target_nii;target];
    CPG1_nii = [CPG1_nii,CPG1];
    CPG2_nii = [CPG2_nii,CPG2];
end


disp('Load data - done')
grad = gradient(sub_pix1(102:1500));
%% Calculate transformation
S4 = [sub_pix2(101:100+N) sub_pix1(101:100+N) ones(N,1)];
S5 = [sub_pix2(101:100+N).^2 sub_pix1(101:100+N) ones(N,1)];
S6 = [grad(101:100+N) sub_pix1(101:100+N) ones(N,1)];
S7 = [grad(101:100+N).^2 sub_pix1(101:100+N) ones(N,1)];
x_1st = S4*c_lin4;
x_2nd = S5*c_2nd5;
x_3rd = S6*c_lin6;
x_4rd = S6*c_2nd7;

disp('Calculate transformation - done')
%% save into registrations
[SI_reg1_1st,AP_reg1_1st,SI_reg2_1st,AP_reg2_1st] = deal(zeros(N,N_CP));
[SI_reg1_2nd,AP_reg1_2nd,SI_reg2_2nd,AP_reg2_2nd] = deal(zeros(N,N_CP));
[SI_reg1_3rd,AP_reg1_3rd,SI_reg2_3rd,AP_reg2_3rd] = deal(zeros(N,N_CP));
[SI_reg1_4rd,AP_reg1_4rd,SI_reg2_4rd,AP_reg2_4rd] = deal(zeros(N,N_CP));

CPG1_nii_1st = CPG1_nii;
CPG2_nii_1st = CPG2_nii;
CPG1_nii_2nd = CPG1_nii;
CPG2_nii_2nd = CPG2_nii;
CPG1_nii_3rd = CPG1_nii;
CPG2_nii_3rd = CPG2_nii;
CPG1_nii_4rd = CPG1_nii;
CPG2_nii_4rd = CPG2_nii;

% X = [SI_reg1_def AP_reg1_def SI_reg2_def AP_reg2_def];
for i=1:N
    
    % linear model
    SI_reg1_1st(i,:) = x_1st(i,(1:N_CP));
    AP_reg1_1st(i,:) = x_1st(i,(1+N_CP:2*N_CP));
    SI_reg2_1st(i,:) = x_1st(i,(1+2*N_CP:3*N_CP));
    AP_reg2_1st(i,:) = x_1st(i,(1+3*N_CP:4*N_CP));
    CPG1_nii_1st(i).img(:,:,1,1,2) = reshape(SI_reg1_1st(i,:),[67,67]);
    CPG1_nii_1st(i).img(:,:,1,1,1) = reshape(AP_reg1_1st(i,:),[67,67]);
    CPG2_nii_1st(i).img(:,:,1,1,2) = reshape(SI_reg2_1st(i,:),[67,67]);
    CPG2_nii_1st(i).img(:,:,1,1,1) = reshape(AP_reg2_1st(i,:),[67,67]);
    
    % second order model
    SI_reg1_2nd(i,:) = x_2nd(i,(1:N_CP));
    AP_reg1_2nd(i,:) = x_2nd(i,(1+N_CP:2*N_CP)); 
    SI_reg2_2nd(i,:) = x_2nd(i,(1+2*N_CP:3*N_CP)); 
    AP_reg2_2nd(i,:) = x_2nd(i,(1+3*N_CP:4*N_CP));
    CPG1_nii_2nd(i).img(:,:,1,1,2) = reshape(SI_reg1_2nd(i,:),[67,67]);
    CPG1_nii_2nd(i).img(:,:,1,1,1) = reshape(AP_reg1_2nd(i,:),[67,67]);
    CPG2_nii_2nd(i).img(:,:,1,1,2) = reshape(SI_reg2_2nd(i,:),[67,67]);
    CPG2_nii_2nd(i).img(:,:,1,1,1) = reshape(AP_reg2_2nd(i,:),[67,67]);
    
    % third order model
    SI_reg1_3rd(i,:) = x_3rd(i,(1:N_CP));
    AP_reg1_3rd(i,:) = x_3rd(i,(1+N_CP:2*N_CP));  
    SI_reg2_3rd(i,:) = x_3rd(i,(1+2*N_CP:3*N_CP));  
    AP_reg2_3rd(i,:) = x_3rd(i,(1+3*N_CP:4*N_CP));
    
    CPG1_nii_3rd(i).img(:,:,1,1,2) = reshape(SI_reg1_3rd(i,:),[67,67]);
    CPG1_nii_3rd(i).img(:,:,1,1,1) = reshape(AP_reg1_3rd(i,:),[67,67]);
    CPG2_nii_3rd(i).img(:,:,1,1,2) = reshape(SI_reg2_3rd(i,:),[67,67]);
    CPG2_nii_3rd(i).img(:,:,1,1,1) = reshape(AP_reg2_3rd(i,:),[67,67]);
    
    SI_reg1_4rd(i,:) = x_4rd(i,(1:N_CP));
    AP_reg1_4rd(i,:) = x_4rd(i,(1+N_CP:2*N_CP));  
    SI_reg2_4rd(i,:) = x_4rd(i,(1+2*N_CP:3*N_CP));  
    AP_reg2_4rd(i,:) = x_4rd(i,(1+3*N_CP:4*N_CP));
    
    CPG1_nii_4rd(i).img(:,:,1,1,2) = reshape(SI_reg1_4rd(i,:),[67,67]);
    CPG1_nii_4rd(i).img(:,:,1,1,1) = reshape(AP_reg1_4rd(i,:),[67,67]);
    CPG2_nii_4rd(i).img(:,:,1,1,2) = reshape(SI_reg2_4rd(i,:),[67,67]);
    CPG2_nii_4rd(i).img(:,:,1,1,1) = reshape(AP_reg2_4rd(i,:),[67,67]);
end

disp('Save into registrations - done')
%% present results & Residual motion
F(N) = struct('cdata',[],'colormap',[]);
Res_motion4 = zeros(N,160,160);
Res_motion5 = zeros(N,160,160);
Res_motion6 = zeros(N,160,160);
Res_motion7 = zeros(N,160,160);

Def_error4=zeros(160,160,1,1,2,N);
Def_error5=zeros(160,160,1,1,2,N);
Def_error6=zeros(160,160,1,1,2,N);
Def_error7=zeros(160,160,1,1,2,N);

for i=1:N
    % plot linear model
    [def_vol_nii_1, def_field_nii_1, dis_field_nii_1] = deformNiiWithCPGsSliding(CPG1_nii_1st(i),CPG2_nii_1st(i),...
        dist_nii,source_nii,target_nii(i),false);
    subplot(2,2,1)
    dispNiiSlice(def_vol_nii_1,'z',1,[],[],[],[],true,[]);
    title('Model 4')
    
    % plot 2nd order model
    [def_vol_nii_2, def_field_nii_2, dis_field_nii_2] = deformNiiWithCPGsSliding(CPG1_nii_2nd(i),CPG2_nii_2nd(i),...
        dist_nii,source_nii,target_nii(i),false);
    subplot(2,2,2)
    dispNiiSlice(def_vol_nii_2,'z',1,[],[],[],[],true,[]);
    title('Model 5')
    
    % plot 3rd order model
    [def_vol_nii_3, def_field_nii_3, dis_field_nii_3] = deformNiiWithCPGsSliding(CPG1_nii_3rd(i),CPG2_nii_3rd(i),...
        dist_nii,source_nii,target_nii(i),false);
    subplot(2,2,3)
    dispNiiSlice(def_vol_nii_3,'z',1,[],[],[],[],true,[]);
    title('Model 6')
    
    [def_vol_nii_4, def_field_nii_4, dis_field_nii_4] = deformNiiWithCPGsSliding(CPG1_nii_4rd(i),CPG2_nii_4rd(i),...
        dist_nii,source_nii,target_nii(i),false);
    subplot(2,2,4)
    dispNiiSlice(def_vol_nii_4,'z',1,[],[],[],[],true,[]);
    title('Model 7')
    
    % save frames
    F(i) = getframe(gcf);
    
    % Residual motion
    Res_motion4(i,:,:) = def_vol_nii_1.img-def_vol_nii.img;
    Res_motion5(i,:,:) = def_vol_nii_2.img-def_vol_nii.img;
    Res_motion6(i,:,:) = def_vol_nii_3.img-def_vol_nii.img;
    Res_motion7(i,:,:) = def_vol_nii_4.img-def_vol_nii.img;

    Def_error4(:,:,1,1,:,i)=def_field_nii_1.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
    Def_error5(:,:,1,1,:,i)=def_field_nii_2.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
    Def_error6(:,:,1,1,:,i)=def_field_nii_3.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
    Def_error7(:,:,1,1,:,i)=def_field_nii_4.img(:,:,1,1,:)-def_field_nii.img(:,:,1,1,:);
    
    
end
SaveAsVideo(F,'myVideo.mp4');


%% Residual motion Results
mu_res_4 = zeros(N,1);
mu_res_5 = zeros(N,1);
mu_res_6 = zeros(N,1);
mu_res_7 = zeros(N,1);
for i=1:N
     res_img_1 = reshape(Res_motion4(i,:,:),[160,160]);
     res_img_2 = reshape(Res_motion5(i,:,:),[160,160]);
     res_img_3 = reshape(Res_motion6(i,:,:),[160,160]);
     res_img_4 = reshape(Res_motion7(i,:,:),[160,160]);
    subplot(2,2,1)
    imshow(res_img_1');colormap('gray');title('Model 4');
    subplot(2,2,2)
    imshow(res_img_2');colormap('gray');title('Model 5');
    subplot(2,2,3)
    imshow(res_img_3');colormap('gray');title('Model 6');
    subplot(2,2,4)
    imshow(res_img_4');colormap('gray');title('Model 7');
    mu_res_4(i) = nansum(nansum(abs(res_img_1)))/(160^2);
    mu_res_5(i) = nansum(nansum(abs(res_img_2)))/(160^2);
    mu_res_6(i) = nansum(nansum(abs(res_img_3)))/(160^2);
    mu_res_7(i) = nansum(nansum(abs(res_img_4)))/(160^2);
end
% figure
% hold on
% plot(mu_res_1,'r-');
% plot(mu_res_2,'g-');
% plot(mu_res_3,'b-');
% plot(mu_res_4,'c-');
% plot(mu_res_5,'m-');
% plot(mu_res_6,'y-');
% plot(mu_res_7,'k-');
% legend('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7')
% title('Residual Motion')
% % saveas(gcf,'plots of residual motion.png')

save('Def_error4','Def_error4');
save('Def_error5','Def_error5');
save('Def_error6','Def_error6');
save('Def_error7','Def_error7');

disp('Plots & Residual motion - done');

disp('####### Question 6.5 task 1 completed #######');
