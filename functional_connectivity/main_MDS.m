clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI  = 'G:\study3_CFS_fMRI_v2\data\mirror_loc_v2.voi';
INPUTDIR       = 'G:\study3_CFS_fMRI_v2\data\';
INPUTPAT_VTC   = '*_msk*mm.vtc'; % only msk vtcs (*msk*mm.vtc)
INPUTPAT_PRT   = '*_msk*vol.prt';
conds          = [2,3]; % msk - 2,3 nmsk - 4,5
tmp_from_onset = 2;
VOIs2Clust     = 1:9;  % excluding outliers based on PDIST results 


voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC, struct('maxdepth',3));
prtfiles = findFilesBVQX(INPUTDIR, INPUTPAT_PRT, struct('maxdepth',4));

% LOAD DATA:
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    prt = BVQXfile(prtfiles{vt});
    for vo = VOIs2Clust % 1:size(voifile.VOI,2)
        if vt==1
            voi_list{vo} = voifile.VOI(vo).Name;
        end
        voi_coords = tal2bv(voifile.VOI(vo).Voxels)';
        voi_rtc(:,vo) = zscore(vtc.VOITimeCourseOrig(voi_coords));
        % cut and concat rtcs:
        low_tmps  = prt.Cond(conds(1)).OnOffsets+tmp_from_onset;
        high_tmps = prt.Cond(conds(2)).OnOffsets+tmp_from_onset;
        voi_rtc_low(:,vo)  = voi_rtc(low_tmps,vo);
        voi_rtc_high(:,vo) = voi_rtc(high_tmps,vo);   
    end
    
    % for each subject:
    rdm(vt,:,:) = squareform(pdist(voi_rtc','correlation')); % corrcoef distance matrix
    rdm_low(vt,:,:)  = squareform(pdist(voi_rtc_low','correlation'));  % corrcoef distance matrix
    rdm_high(vt,:,:) = squareform(pdist(voi_rtc_high','correlation')); % corrcoef distance matrix
    
    clear voi_rtc
    vtc.ClearObject();
end
mrdm      = mean(rdm);
mrdm      = reshape(mrdm, size(mrdm,2), size(mrdm,3));
mrdm_low  = mean(rdm_low);
mrdm_high = mean(rdm_high);
mrdm_low  = reshape(mrdm_low, size(mrdm_low,2), size(mrdm_low,3));
mrdm_high = reshape(mrdm_high, size(mrdm_high,2), size(mrdm_high,3));


%% 1. PLOT 2D MDS (look at the distances):
[Y,e] = cmdscale(1-mrdm_high);
plot(-Y(:,1),Y(:,2),'.');
text(-Y(:,1),Y(:,2),voi_list);



%% 3. single ROIs STATISTICS
% for s = unique(sub_vtcs_idx)
%     for r = 1:12
%         rois_mean_dist(s,r) = sum(sub_rdm(s,r,:))/12;
%     end
% end

%convert 9*9 rdm to 8*8 rdm (average PMds)
mrdm2convert = mrdm_high;
pmds = mrdm2convert(:,[3 4]);
mean_pmd9 = mean(pmds,2);
mean_pmd_vertical = mean(mean_pmd9([3,4]));
mean_pmd8 = [mean_pmd9(1:2); mean_pmd_vertical; mean_pmd9(5:end)];
new_mat98 = [mrdm2convert(:,[1 2]), mean_pmd9, mrdm2convert(:,5:end)]; 
new_mat88 = [new_mat98([1 2],:); mean_pmd8'; new_mat98(5:end,:)];
new_mat88(3,3)=0;
% figure; imagesc(1-new_mat88, [0 0.6])


    

