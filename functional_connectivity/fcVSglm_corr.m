clear all; close all; clc;
INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR_VTC  = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_msk*mm.vtc'; % only msk vtcs

voi_file  = BVQXfile(INPUTFILE_VOI);
vtc_files = findFilesBVQX(INPUTDIR_VTC,INPUTPAT_VTC,struct('maxdepth',3));

%% compute FC matrix
for vt = 1:length(vtc_files)
    fprintf('loading vtc %1.0f / %2.0f\n', vt, length(vtc_files))
    vtc = BVQXfile(vtc_files{vt});
    for vo1 = 1:12; % length(voi_file.VOI)
        voi_list{vo1} = voi_file.VOI(vo1).Name;
        % fprintf('processing voi %1.0f / %2.0f\n', vo1, length(voi_file.VOI))
        voi1_coords = tal2bv(voi_file.VOI(vo1).Voxels)';
        voi1_rtc = zscore(vtc.VOITimeCourseOrig(voi1_coords));
        for vo2 = 1:12; %length(voi_file.VOI)
            voi2_coords = tal2bv(voi_file.VOI(vo2).Voxels)';
            voi2_rtc = zscore(vtc.VOITimeCourseOrig(voi2_coords));
            [r,p] = corrcoef(voi1_rtc,voi2_rtc);
            r_mat(vo1,vo2,vt) = r(2);
        end
    end
    vtc.ClearObject(); clear vtc
end
meanr_mat = mean(r_mat,3);
meanr_mat(logical(eye(size(meanr_mat)))) = 0;


%% comput GLM matrix
load('glm_M-NM_vals.mat');
for vo1 = 1:12; % length(voi_file.VOI)
    % fprintf('processing voi %1.0f / %2.0f\n', vo1, length(voi_file.VOI))
    voi1_val = vals{vo1,2}(1);
    for vo2 = 1:12; %length(voi_file.VOI)
        voi2_val = vals{vo2,2}(1);
        meand = abs(voi1_val - voi2_val); % M-NM diff matrix
        meand_mat(vo1,vo2) = meand;
    end
end
r_arr = squareform(meanr_mat);
d_arr = squareform(meand_mat);

[r,p] = corrcoef(d_arr, r_arr);
fit = polyfit(d_arr, r_arr, 1);
line = fit(1)*d_arr + fit(2)*ones(1,size(d_arr,2));

scatter(d_arr, r_arr); hold on;
plot(d_arr, line, 'b'); hold on;
xlabel('ROIs beta difference (M(H-L) - NM(H-L))');
ylabel('ROIs correlation coeffitient (R)');
title(['without visual ROIs: r = ' num2str(r(1,2)) ', pval = ' num2str(p(1,2))]);




