%% Q6_5_2 Deformation Field Error
%Attention:firslty run Q6_5_1, then this!

%The deformation field is calculated by euclidean distance
%of AP(which could be regarded as delta(x)) and SI(which could be regarded as delta(y))

% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

N=50;%number of t-frame set which must be the same N as in Q6.5.1
N_r=160;
N_c=160;
N_pixel = N_r*N_c;

load('Def_error4.mat');
load('Def_error5.mat');
load('Def_error6.mat');
load('Def_error7.mat');

%Record the L2 of AP deformation error, SI deformation error and the total
%deformation error over time for every pixel.
[L2_AP_1st,L2_AP_2nd,L2_AP_3rd,L2_AP_4rd...
 L2_SI_1st,L2_SI_2nd,L2_SI_3rd,L2_SI_4rd...
 L2_1st,L2_2nd,L2_3rd,L2_4rd]=deal(zeros(N_r,N_c));

%Record the averaged AP deformation error, averaged SI deformation error and the averaged total deformation error of
%over every pixel for every time frame.
[error_AP_1st,error_AP_2nd,error_AP_3rd,error_AP_4rd...
 error_SI_1st,error_SI_2nd,error_SI_3rd,error_SI_4rd...
 error_4st,error_5nd,error_6rd,error_7rd]=deal(zeros(N,1));

source_mask = load_untouch_nii('0007_mask.nii');
dist_nii = load_untouch_nii('0007_sdt.nii');
target_mask=zeros(N_r,N_c,N);

%Before evaluation in the target image space, the mask images in source
%space needed to be transformed by the transformation found by the
%registrations.
for i=101:100+N 
    if i<10
        f_name = sprintf('000%i', i);
        elseif (10<=i) && (i<100)
        f_name = sprintf('00%i', i);
        elseif  (100<=i) && (i<1000)
        f_name = sprintf('0%i', i);
    else
        f_name = sprintf('%i', i);
    end
    
    %original images
    nii = load_untouch_nii(sprintf('%s.nii', f_name));
    target_nii=nii;
    %The motion of region 1
    CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    %The motion of region 2
    CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
    %To transform the 0007mask from source image space to target image space by
    %the transformation found by the registration.
    [target_struct, ~,~] = deformNiiWithCPGsSliding(CPG1_nii,CPG2_nii,...
    dist_nii,source_mask,target_nii,false);
    target_mask(:,:,i-100)=target_struct.img;
    
    %eliminate the NaN and numbers which are not 0 and 1
    for r=1:160
        for c=1:160
            if target_mask(r,c,i-100)>0.5
               target_mask(r,c,i-100)=1;
            else 
               target_mask(r,c,i-100)=0;
            end 
        end 
    end 
  fprintf('target space mask information of %i th frame has been stored.\n',i);  
 end 
%calculate L2 norm of AP deformation field error,SI deformation field error at each pixel of N time series for 3 models
%respectively
for r=1:N_r
    for c=1:N_c
           n=0;
           mask_this_pixel=reshape(target_mask(r,c,:),1,[]);
           
           L2_AP_1st(r,c)=nansum((reshape(Def_error4(r,c,1,1,1,:),1,[]).^2).*mask_this_pixel);
           L2_AP_2nd(r,c)=nansum((reshape(Def_error5(r,c,1,1,1,:),1,[]).^2).*mask_this_pixel);
           L2_AP_3rd(r,c)=nansum((reshape(Def_error6(r,c,1,1,1,:),1,[]).^2).*mask_this_pixel);
           L2_AP_4rd(r,c)=nansum((reshape(Def_error7(r,c,1,1,1,:),1,[]).^2).*mask_this_pixel);
           
           L2_SI_1st(r,c)=nansum((reshape(Def_error4(r,c,1,1,2,:),1,[]).^2).*mask_this_pixel);
           L2_SI_2nd(r,c)=nansum((reshape(Def_error5(r,c,1,1,2,:),1,[]).^2).*mask_this_pixel);
           L2_SI_3rd(r,c)=nansum((reshape(Def_error6(r,c,1,1,2,:),1,[]).^2).*mask_this_pixel);
           L2_SI_4rd(r,c)=nansum((reshape(Def_error7(r,c,1,1,2,:),1,[]).^2).*mask_this_pixel);
           
           L2_1st(r,c)= L2_AP_1st(r,c)+ L2_SI_1st(r,c);
           L2_2nd(r,c)= L2_AP_2nd(r,c)+ L2_SI_2nd(r,c);
           L2_3rd(r,c)= L2_AP_3rd(r,c)+ L2_SI_3rd(r,c);
           L2_4rd(r,c)= L2_AP_4rd(r,c)+ L2_SI_4rd(r,c);
           fprintf('The L2 norm of  AP,SI error of pixel(%i,%i) has been obtained.\n',r,c);
    end
    
end 
%calculate  AP deformation field error,SI deformation field error at each frame(totally N time series) for 3 models
%respectively and combined deformation deformation field error by euclidean distance at each frame(totally N time series) for 3 models
%respectively
for t=1:N
           mask_this_frame=reshape(target_mask(:,:,t),1,[]);
           
           error_AP_1st(t)=nansum(reshape(abs(Def_error4(:,:,1,1,1,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_AP_2nd(t)=nansum(reshape(abs(Def_error5(:,:,1,1,1,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_AP_3rd(t)=nansum(reshape(abs(Def_error6(:,:,1,1,1,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_AP_4rd(t)=nansum(reshape(abs(Def_error7(:,:,1,1,1,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           
           error_SI_1st(t)=nansum(reshape(abs(Def_error4(:,:,1,1,2,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_SI_2nd(t)=nansum(reshape(abs(Def_error5(:,:,1,1,2,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_SI_3rd(t)=nansum(reshape(abs(Def_error6(:,:,1,1,2,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           error_SI_4rd(t)=nansum(reshape(abs(Def_error7(:,:,1,1,2,t)),1,[]).*mask_this_frame)/(N_r*N_c);
           
           error_4st(t)=sqrt(error_AP_1st(t)^2+error_SI_1st(t)^2);
           error_5nd(t)=sqrt(error_AP_2nd(t)^2+error_SI_2nd(t)^2);
           error_6rd(t)=sqrt(error_AP_3rd(t)^2+error_SI_3rd(t)^2);
           error_7rd(t)=sqrt(error_AP_4rd(t)^2+error_SI_4rd(t)^2);

           fprintf('The AP,SI error of time %i has been obtained.\n',t);
end
% error_AP_1st=sum(error_AP_1st)/N
% error_AP_2nd=sum(error_AP_2nd)/N
% error_AP_3rd=sum(error_AP_3rd)/N
% 
% error_SI_1st=sum(error_SI_1st)/N
% error_SI_2nd=sum(error_SI_2nd)/N
% error_SI_3rd=sum(error_SI_3rd)/N
% 
% error_1st=error_1st;
% error_2nd=sum(error_2nd);
% error_3rd=sum(error_3rd);
figure;
%show the result of the deformation field error (Q6_5_2)
plot(error_1st,'r-');
hold on;
plot(error_2nd,'g-');
hold on;
plot(error_3rd,'b-');
hold on;
plot(error_4st,'c-');
hold on;
plot(error_5nd,'m-');
hold on;
plot(error_6rd,'y-');
hold on;
plot(error_7rd,'k-');
legend('Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7')
title('Deformation Field Errors')
% %show the result of the visual assessment (Q6_5_1)
% subplot(1,2,2);
% plot(mu_res_1);
% hold on;
% plot(mu_res_2);
% hold on;
% plot(mu_res_3);
% legend('linear','2nd poly','3rd poly')
% title('Q6_5_1:residual motion')

%The L2 norm of the deformation field error at each pixel(AP and SI
%seperately) inside of the deformation field masks.
save('L2_AP_1st','L2_AP_1st');
save('L2_AP_2nd','L2_AP_2nd');
save('L2_AP_3rd','L2_AP_3rd');

save('L2_SI_1st','L2_SI_1st');
save('L2_SI_2nd','L2_SI_2nd');
save('L2_SI_3rd','L2_SI_3rd');

save('L2_1st','L2_1st');
save('L2_2nd','L2_2nd');
save('L2_3rd','L2_3rd');

