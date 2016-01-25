%% ppi_results_plotsAndStat

load('ppiMasked_VOIGLM_manual_voisorted.mat');    ma_mat = beta_mat;
load('ppiNonMasked_VOIGLM-manual_voisorted.mat'); nm_mat = beta_mat;

delta_ma = ma_mat(:,:,:,2) - ma_mat(:,:,:,1); % masked High > Low
delta_nm = nm_mat(:,:,:,2) - nm_mat(:,:,:,1); % non-masked High > Low

effect_size_mat = delta_ma - delta_nm;
avg_effect_size = mean(effect_size_mat(:,:,:,1),1);
avg_effect_size = reshape(avg_effect_size, size(avg_effect_size,2), ...
    size(avg_effect_size,3));

% Statistics:
for vo1 = 1:size(delta_ma,2)
    for vo2 = 1:size(delta_ma,3)
        [~,p(vo1,vo2)] = ttest(delta_ma(:,vo1,vo2), delta_nm(:,vo1,vo2), 'tail', 'right');
    end
end
[~, crit_p, ~, adj_p]=fdr_bh(p(:));
q = reshape(adj_p, size(p,1), size(p,2));
q(q > 0.1) = 1;

subplot(1,2,1)
imagesc(avg_effect_size);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voi_list));
set(gca,'YTick',1:length(voi_list));
set(gca,'XAxisLocation','top');
title('delta M(H>L) > delta NM(H>L) BETAS');
caxis([-0.1, 0.1]);

subplot(1,2,2)
imagesc(q);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voi_list));
set(gca,'YTick',1:length(voi_list));
set(gca,'XAxisLocation','top');
title('delta M(H>L) > delta NM(H>L) Pvalues');
caxis([-0.2, 0.2]);

% TODO: 
% ttests for each voi*voi + fdr correction