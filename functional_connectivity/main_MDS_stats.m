%% 2. STATISTICS
% single subjects rdms
load('sub_vtcs_idx.mat');
load('rdm_rest.mat');
sub_rdm = sub_rdm_rest;
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