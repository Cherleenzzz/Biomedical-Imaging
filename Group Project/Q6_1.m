% add path
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/NIfTI_20140122');
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/display')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/images')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/registrations')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/CMBI-project/transforms')
addpath('/Users/humingzhou/Downloads/UCL/Biomedical_Image/Project/segmentation')

%% images
% load and view .nii files
nii=cell(1500,1);
figure
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
    nii{i} = load_untouch_nii(sprintf('%s.nii', f_name));
    dispNiiSlice(nii{i},'z',1,[],[],[],[],true,[])
    title('MR image')
    pause(0.001)
end

%% registration
% load registration results
dist_nii = load_untouch_nii('0007_sdt.nii');
source_nii = load_untouch_nii('0007.nii');
figure
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
    CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
    CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
    [def_vol_nii, def_field_nii, dis_field_nii] = deformNiiWithCPGsSliding(CPG1_nii,CPG2_nii,...
        dist_nii,source_nii,nii{i},false);
    dispNiiSlice(def_vol_nii,'z',1,[],[],[],[],true,[])
    title('Registration result')
    pause(0.001)
end
% 
% %% display results side by side
% figure
% for i=1:1500
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
%     % original images
%     nii = load_untouch_nii(sprintf('%s.nii', f_name));
%     subplot(1,2,1)
%     dispNiiSlice(nii,'z',1,[],[],[],[],true,[])
%     title('original')
%     
%     % registration results
%     CPG1_nii = load_untouch_nii(sprintf('%s_cpp_region1.nii', f_name));
%     CPG2_nii = load_untouch_nii(sprintf('%s_cpp_region2.nii', f_name));
%     
%     [def_vol_nii, def_field_nii, dis_field_nii] = deformNiiWithCPGsSliding(CPG1_nii,CPG2_nii,...
%         dist_nii,source_nii,nii{i},false);
%     subplot(1,2,2)
%     dispNiiSlice(def_vol_nii,'z',1,[],[],[],[],true,[])
%     title('registered')
%     pause(0.001)
% end

