
%% TODO
%  p val across subjects
%  k means on voi rtcs with k=2
%  rsa? 

clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final_ppisorted.voi';
INPUTDIR      = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_msk*mm.vtc'; % only msk vtcs
INPUTPAT_SDM  = '*_msk*mm*.sdm';
OUTPUTDIR_DMS = 'D:\study 3_CFS-fMRI_v2\analysis\functional_connectivity\ppi_dms';
PRD_LIST      = '"seedvoi*ML" "seedvoi*MH" "prdML" "prdMH" "seedvoi" "constant"';
GLM_NAME      = 'manual_voisorted';
GENERATE_DMS  = 0; % generate sdms and mdms
GENERATE_GLM  = 1;
CONTRASTS     = [1 0 0 0 0 0; 0 1 0 0 0 0]; % contrasts for not-confound predictors

prdAcol       = 1; % colunm number in sdm file of condition A
prdBcol       = 2; % colunm number in sdm file of condition B
constant_col  = 7; % colunm number in sdm file of the constant

voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC,struct('maxdepth',3));
sdmfiles = findFilesBVQX(INPUTDIR, INPUTPAT_SDM,struct('maxdepth',3));

%% MAIN
if GENERATE_GLM
    for vo1 = 1:size(voifile.VOI,2)
        %% 1. get seed voi rtc
        seed_coords = tal2bv(voifile.VOI(vo1).Voxels)';
        seed_name   = voifile.VOI(vo1).Name;
        fprintf('\n ppi processing for seed:  %s' , seed_name);
        
        for vt = 1:length(vtcfiles)
            vtc = BVQXfile(vtcfiles{vt});
            vtcnamepart = strsplit(vtcfiles{vt}, '\');
            vtcnamepart = strsplit(vtcnamepart{6}, '.');
            vtcnamepart = vtcnamepart{1};
            seed_rtc = zscore(vtc.VOITimeCourseOrig(seed_coords));
            
            if GENERATE_DMS
                %% 1. define interaction predictors
                seedA = zscore(seed_rtc) .* prdA ;
                seedB = zscore(seed_rtc) .* prdB ;
                all_prds = [seedA, seedB, prdA, prdB, seed_rtc, constant];
                
                %% 2. define main predictors
                sdm = BVQXfile(sdmfiles{vt});
                prdA = zscore(sdm.SDMMatrix(:,prdAcol));
                prdB = zscore(sdm.SDMMatrix(:,prdBcol));
                constant = sdm.SDMMatrix(:,constant_col);
                
                %% 3. write sdm file
                fid = fopen([OUTPUTDIR_DMS '\' seed_name '_'  vtcnamepart '.sdm'],'wt');
                % sdm header :
                fprintf(fid,'\n%s\n','FileVersion:            1');
                fprintf(fid,'\n%s\n','NrOfPredictors:         6');
                fprintf(fid,'\n%s%\n',['NrOfDataPoints:         ', num2str(size(vtc.VTCData,1))]);
                fprintf(fid,'\n%s\n','IncludesConstant:       0') ;
                fprintf(fid,'\n%s\n','FirstConfoundPredictor: 3') ;
                fprintf(fid,'\n%s\n','255 0 255   0 255 255   255 0 0   0 255 0   0 0 255') ;
                fprintf(fid,'\n%s\n', PRD_LIST) ;
                for row = 1:size(all_prds,1)
                    string = sprintf('%6.6f\t',all_prds(row,:));
                    fprintf(fid,'%s\r\n',string);
                end
                fclose(fid);
                sdm.ClearObject; clear sdm;
            end
            
            %% 5. do manual VOIGLM and extract voi's betas
            if GENERATE_GLM
                ppi_sdm = BVQXfile([OUTPUTDIR_DMS '\' seed_name '_'  vtcnamepart '.sdm']);
                [b, p]  = voiglm(ppi_sdm, vtc, CONTRASTS', voifile);
                beta_mat(vt,vo1,:,:) = b; % (vtc, voi1(seed), voi2, contrast)
                pval_mat(vt,vo1,:,:) = p;
                vtc.ClearObject; clear vtc;
            end
        end
        voi_list{vo1} = seed_name;
    end
    save(['ppi_VOIGLM-' GLM_NAME '.mat'], 'beta_mat', 'pval_mat', 'CONTRASTS', 'voi_list');
    
else
    load(['ppi_VOIGLM-' GLM_NAME   '.mat'])
end


%% plot results
meanb_ML_mat = mean(beta_mat(:,:,:,1),1);
meanb_MH_mat = mean(beta_mat(:,:,:,2),1);
meanp_ML_mat = mean(pval_mat(:,:,:,1),1);
meanp_MH_mat = mean(pval_mat(:,:,:,2),1);
meanp_ML_mat(meanp_ML_mat > 0.1) = 1;
meanp_MH_mat(meanp_MH_mat > 0.1) = 1;

meanb_ML_mat = reshape(meanb_ML_mat, size(meanb_ML_mat,2), size(meanb_ML_mat,3));
meanb_MH_mat = reshape(meanb_MH_mat, size(meanb_MH_mat,2), size(meanb_MH_mat,3));
meanp_ML_mat = reshape(meanp_ML_mat, size(meanp_ML_mat,2), size(meanp_ML_mat,3));
meanp_MH_mat = reshape(meanp_MH_mat, size(meanp_MH_mat,2), size(meanp_MH_mat,3));

subplot(2,2,1)
imagesc(meanb_MH_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('Mean Beta - MSK HIGH');
caxis([-0.0675, 0.0675]); 

subplot(2,2,2)
imagesc(meanp_MH_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('P Val - MSK HIGH');

subplot(2,2,3)
imagesc(meanb_ML_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('Mean Beta - MSK LOW');
caxis([-0.0675, 0.0675]); 

subplot(2,2,4)
imagesc(meanp_ML_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('P Val - MSK LOW');








