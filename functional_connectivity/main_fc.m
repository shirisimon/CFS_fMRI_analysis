%% main
clear all; close all; clc;

INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR_VTC  = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_msk*mm.vtc'; % only msk vtcs

voi_file  = BVQXfile(INPUTFILE_VOI);
vtc_files = findFilesBVQX(INPUTDIR_VTC,INPUTPAT_VTC,struct('maxdepth',3));

for vt = 1:length(vtc_files)
    fprintf('loading vtc %1.0f / %2.0f\n', vt, length(vtc_files))
    vtc = BVQXfile(vtc_files{vt});
    for vo1 = 1:length(voi_file.VOI)
        voi_list{vo1} = voi_file.VOI(vo1).Name;
        % fprintf('processing voi %1.0f / %2.0f\n', vo1, length(voi_file.VOI))
        voi1_coords = tal2bv(voi_file.VOI(vo1).Voxels)';
        voi1_rtc = zscore(vtc.VOITimeCourseOrig(voi1_coords));
        for vo2 = 1:length(voi_file.VOI)
            voi2_coords = tal2bv(voi_file.VOI(vo2).Voxels)';
            voi2_rtc = zscore(vtc.VOITimeCourseOrig(voi2_coords));
            [r,p] = corrcoef(voi1_rtc,voi2_rtc);
            r_mat(vo1,vo2,vt) = r(2);
            p_mat(vo1,vo2,vt) = p(2);
        end
    end
    vtc.ClearObject(); clear vtc
end



%% plot results
meanr_mat = mean(r_mat,3);
meanp_mat = mean(p_mat,3);
meanp_mat(meanp_mat > 0.05) = 1;

subplot(1,2,1)
imagesc(meanr_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7.5)
set(gca,'XTick',1:length(voi_file.VOI));
set(gca,'YTick',1:length(voi_file.VOI));
set(gca,'XAxisLocation','top');
title('Corr Coef (pearson r)');

subplot(1,2,2)
imagesc(meanp_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7.5)
set(gca,'XTick',1:length(voi_file.VOI));
set(gca,'YTick',1:length(voi_file.VOI));
set(gca,'XAxisLocation','top');
title('Corr Sig (p val)');

