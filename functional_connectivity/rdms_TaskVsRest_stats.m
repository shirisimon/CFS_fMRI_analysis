%% rdms_stats

rdm_task = importdata('rdm_allconds.mat');
sub_rdm_rest = importdata('rdm_rest.mat');

% downsample rdm_rest to subjects:
load('sub_vtcs_idx.mat');
for s = unique(sub_vtcs_idx)
    [~,sub_idx] = find(sub_vtcs_idx==s);
    sub_rdm_task(s,:,:) = mean(rdm_task(sub_idx,:,:));
end

for i=1:size(sub_rdm_rest,2)
    for j = 1:size(sub_rdm_rest,3)
        [~,pval(i,j)] = ttest2(sub_rdm_rest(:,i,j), sub_rdm_task(:,i,j), ...
            'tail', 'both', 'vartype', 'unequal');
        if i==j
            pval(i,j) = 0;
        end
    end
end

% correct sig_rdm to fdr
pval_vec = squareform(pval);
