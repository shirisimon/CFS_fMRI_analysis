%% ppi_results_plotsAndStat

load('ppi_VOIGLM-MA_manual_voisorted.mat'); ma_mat = beta_mat;
load('ppi_VOIGLM-NM_manual_voisorted.mat'); nm_mat = beta_mat;

delta_ma = ma_mat(:,:,:,2) - ma_mat(:,:,:,1); % masked High > Low
delta_nm = nm_mat(:,:,:,2) - nm_mat(:,:,:,1); % non-masked High > Low

effect_size_mat = delta_ma - delta_nm;
avg_effect_size = mean(effect_size_mat(:,:,:,1),1);
avg_effect_size = reshape(avg_effect_size, size(avg_effect_size,2), ...
    size(avg_effect_size,3));

% M(H>L) Statistics: 
for vo1 = 1:size(ma_mat,2)
    for vo2 = 1:size(ma_mat,3)
        [~,p_ma(vo1,vo2)] = ttest(ma_mat(:,vo1,vo2,2), ma_mat(:,vo1,vo2,1), 'tail', 'right');
    end
end
[~, ~, ~, adj_p]=fdr_bh(p_ma(:));
q_ma = reshape(adj_p, size(p_ma,1), size(p_ma,2));
p_ma(p_ma > 0.05) = 1;
q_ma(q_ma > 0.05) = 1;

% M(H>L) Statistics: 
for vo1 = 1:size(nm_mat,2)
    for vo2 = 1:size(nm_mat,3)
        [~,p_nm(vo1,vo2)] = ttest(nm_mat(:,vo1,vo2,2), nm_mat(:,vo1,vo2,1), 'tail', 'right');
    end
end
[~, ~, ~, adj_p]=fdr_bh(p_nm(:));
q_nm = reshape(adj_p, size(p_nm,1), size(p_nm,2));
p_nm(p_nm > 0.05) = 1;
q_nm(q_nm > 0.05) = 1;


% M>NM Statistics:
for vo1 = 1:size(delta_ma,2)
    for vo2 = 1:size(delta_ma,3)
        [~,p(vo1,vo2)] = ttest(delta_ma(:,vo1,vo2), delta_nm(:,vo1,vo2), 'tail', 'right');
    end
end
[~, crit_p, ~, adj_p]=fdr_bh(p(:));
q = reshape(adj_p, size(p,1), size(p,2));
p(p > 0.05) = 1;
q(q > 0.05) = 1;

delta_ma = mean(delta_ma(:,:,:,1),1);
delta_nm = mean(delta_nm(:,:,:,1),1);
delta_ma = reshape(delta_ma, size(delta_ma,2), size(delta_ma,3));
delta_nm = reshape(delta_nm, size(delta_nm,2), size(delta_nm,3));

 
subplot(1,3,1)
imagesc(delta_ma); % avg_effect_size
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voi_list));
set(gca,'YTick',1:length(voi_list));
set(gca,'XAxisLocation','top');
title('delta M(H>L) BETAS'); % delta M(H>L) > 
caxis([-0.03, 0.03]);

subplot(1,3,2)
imagesc(p_ma);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voi_list));
set(gca,'YTick',1:length(voi_list));
set(gca,'XAxisLocation','top');
title('delta M(H>L) Pvalues');
caxis([-0.03, 0.03]);

subplot(1,3,3)
imagesc(q_ma);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voi_list));
set(gca,'YTick',1:length(voi_list));
set(gca,'XAxisLocation','top');
title('delta M(H>L) Corrected Pvalues');
caxis([-0.03, 0.03]);
