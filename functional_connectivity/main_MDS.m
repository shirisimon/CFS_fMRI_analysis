clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI  = 'G:\study3_CFS_fMRI_v2\data\mirror_loc_v2.voi';
INPUTDIR       = 'G:\study3_CFS_fMRI_v2\data\rest_hcp\';
INPUTPAT_VTC   = '*.vtc'; % only msk vtcs (*msk*mm.vtc)
INPUTPAT_PRT   = '*msk*vol.prt';
conds          = [2,3,4,5]; % msk - 2,3 nmsk - 4,5
tmp_from_onset = 2;
VOIs2Clust     = 1:9;  % excluding outliers based on PDIST results 


voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC, struct('maxdepth',1));
prtfiles = findFilesBVQX(INPUTDIR, INPUTPAT_PRT, struct('maxdepth',4));

% LOAD DATA:
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    %prt = BVQXfile(prtfiles{vt});
    for vo = VOIs2Clust % 1:size(voifile.VOI,2)
        if vt==1
            voi_list{vo} = voifile.VOI(vo).Name;
        end
        voi_coords = tal2bv(voifile.VOI(vo).Voxels)';
        voi_rtc(:,vo) = zscore(vtc.VOITimeCourseOrig(voi_coords));
        % cut and concat rtcs:
        %low_tmps  = prt.Cond(conds(1)).OnOffsets+tmp_from_onset;
        %high_tmps = prt.Cond(conds(2)).OnOffsets+tmp_from_onset;
        %voi_rtc_low(:,vo)  = voi_rtc(low_tmps,vo);
        %voi_rtc_high(:,vo) = voi_rtc(high_tmps,vo);   
    end
    
    % for each subject:
    rdm(vt,:,:) = squareform(pdist(voi_rtc','correlation')); % corrcoef distance matrix
    %rdm_low(vt,:,:)  = squareform(pdist(voi_rtc_low','correlation'));  % corrcoef distance matrix
    %rdm_high(vt,:,:) = squareform(pdist(voi_rtc_high','correlation')); % corrcoef distance matrix
    
    clear voi_rtc
    vtc.ClearObject();
end
mrdm      = mean(rdm);
mrdm      = reshape(mrdm, size(mrdm,2), size(mrdm,3));
% mrdm_low  = mean(rdm_low);
% mrdm_high = mean(rdm_high);
% mrdm_low  = reshape(mrdm_low, size(mrdm_low,2), size(mrdm_low,3));
% mrdm_high = reshape(mrdm_high, size(mrdm_high,2), size(mrdm_high,3));


%% 1. PLOT 2D MDS (look at the distances):
%[Y,e] = cmdscale(mrdm);
%plot(-Y(:,1),Y(:,2),'.');
%text(-Y(:,1),Y(:,2),voi_list);

%% 2. STATISTICS
% single subjects rdms
load('sub_vtcs_idx.mat');
sen_rois = 1:5;
ins_rois = 6:9;
for s = unique(sub_vtcs_idx)
    sub_rdm = rdm;
    %[~,sub_idx] = find(sub_vtcs_idx==s);
    %sub_rdm(s,:,:) = mean(rdm(sub_idx,:,:));
    %sub_rdm_low(s,:,:) = mean(rdm_low(sub_idx,:,:));
    %sub_rdm_high(s,:,:) = mean(rdm_high(sub_idx,:,:));
    
    %% calc mean distance within sensitive regions:
    % all:
    within_sen_mean_dist_all(s,1) = sum(sum(sub_rdm(s,sen_rois,sen_rois)))/...
        length(sen_rois)^2;
    % low:
    %within_sen_mean_dist_low(s,1) = sum(sum(sub_rdm_low(s,sen_rois,sen_rois)))/...
    %    length(sen_rois)^2;
    % high:
    %within_sen_mean_dist_high(s,1) = sum(sum(sub_rdm_high(s,sen_rois,sen_rois)))/...
    %    length(sen_rois)^2;
    
    %% calc mean distance within insensitive regions:
    % all:
    within_ins_mean_dist_all(s,1) = sum(sum(sub_rdm(s,ins_rois,ins_rois)))/...
        length(ins_rois)^2;
    % low:
    %within_ins_mean_dist_low(s,1) = sum(sum(sub_rdm_low(s,ins_rois,ins_rois)))/...
    %    length(ins_rois)^2;
    % high:
    %within_ins_mean_dist_high(s,1) = sum(sum(sub_rdm_high(s,ins_rois,ins_rois)))/...
    %    length(ins_rois)^2;
    
    %% calc mean distance between sensitive and insensitive regions:
    % all
    between_groups_mean_dist_all(s,1) = sum(sum(sub_rdm(s,sen_rois,ins_rois)))/...
        (length(sen_rois)*length(ins_rois));
    % low
    %between_groups_mean_dist_low(s,1) = sum(sum(sub_rdm_low(s,sen_rois,ins_rois)))/...
    %    (length(sen_rois)*length(ins_rois));
    % high
    %between_groups_mean_dist_high(s,1) = sum(sum(sub_rdm_high(s,sen_rois,ins_rois)))/...
    %    (length(sen_rois)*length(ins_rois));
end

%% 3. single ROIs STATISTICS
% for s = unique(sub_vtcs_idx)
%     for r = 1:12
%         rois_mean_dist(s,r) = sum(sub_rdm(s,r,:))/12;
%     end
% end

    

