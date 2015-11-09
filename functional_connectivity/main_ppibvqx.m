%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WARNING: beta values computed in this code are not consistent with BV beta values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR      = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_msk*mm.vtc'; % only msk vtcs
INPUTPAT_SDM  = '*_msk*mm*.sdm';
OUTPUTDIR_DMS = 'D:\study 3_CFS-fMRI_v2\analysis\functional_connectivity\ppi_dms';
PRD_LIST      = '"seedvoi*ML" "seedvoi*MH" "prdML" "prdMH" "seedvoi" "constant"';
GENERATE_DMS  = 1; % generate sdms and mdms
LOAD_MDMS     = 0; % to generate new glm with different parameters (opts)
LOAD_GLM      = 0; % load beta and pval of previous glm results.
GLM_NAME      = 'robust';

prdAcol       = 1; % colunm number in sdm file of condition A
prdBcol       = 2; % colunm number in sdm file of condition B
constant_col  = 7; % colunm number in sdm file of the constant

voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC,struct('maxdepth',3));
sdmfiles = findFilesBVQX(INPUTDIR, INPUTPAT_SDM,struct('maxdepth',3));

%% MAIN
if ~LOAD_GLM
    for vo1 = 1:length(voifile.VOI)
        %% 1. get seed voi rtc
        seed_coords = tal2bv(voifile.VOI(vo1).Voxels)';
        seed_name = voifile.VOI(vo1).Name;
        voi_list{vo1} = seed_name;
        fprintf('\n ppi processing for seed:  %s' , seed_name);
        
        if GENERATE_DMS
            for vt = 1:length(vtcfiles)
                %% 2. define main predictors
                sdm = BVQXfile(sdmfiles{vt});
                prdA = zscore(sdm.SDMMatrix(:,prdAcol));
                prdB = zscore(sdm.SDMMatrix(:,prdBcol));
                constant = sdm.SDMMatrix(:,constant_col);
                
                %% 3. define interaction predictors
                vtc = BVQXfile(vtcfiles{vt});
                seed_rtc = zscore(vtc.VOITimeCourse(seed_coords));
                seedA = zscore(seed_rtc) .* prdA ;
                seedB = zscore(seed_rtc) .* prdB ;
                all_prds = [seedA, seedB, prdA, prdB, seed_rtc, constant] ;
                
                %% 4. write sdm file
                vtcpart = strsplit(vtcfiles{vt}, '\');
                vtcpart = strsplit(vtcpart{6}, '.');
                vtcpart = vtcpart{1};
                fid = fopen([OUTPUTDIR_DMS '\' seed_name '_'  vtcpart '.sdm'],'wt');
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
                vtc.ClearObject; clear vtc;
            end
            ppi_sdmfiles = findFilesBVQX(OUTPUTDIR_DMS, ['*' seed_name '*.sdm']);
            
            %% 4. generate mdm and multi-subject ROI glm
            mdm = xff('new:mdm');
            mdm.RFX_GLM = 1;
            mdm.zTransformation = 1;
            mdm.NrOfStudies = numel(vtcfiles);
            mdm.XTC_RTC = [vtcfiles(:), ppi_sdmfiles(:)];
            mdm.SaveAs([OUTPUTDIR_DMS '\' seed_name '.mdm']);
            
        elseif LOAD_MDMS
            mdmfile = findFilesBVQX(OUTPUTDIR_DMS, [seed_name '.mdm']);
            mdm = BVQXfile(mdmfile{1});
            mdm.RFX_GLM = 1;
        end
        
        %% 5. do VOIGLM and extract voi's betas
         GLM_OPTS = struct('robust', true, ...
                   'ppicond', {'seedvoi*ML', 'seedvoi*ML'}, ...
                   'ppivoi', voifile);
        [~,b] = mdm.ComputeVOIGLM(voifile, GLM_OPTS);
        [~,meanp] = ttest(b);
        beta_mat(:,:,:,vo1) = b;  % (vt, b, vo2, vo1)
        pval_mat(:,:,:,vo1) = meanp;
        
        mdm.ClearObject; clear mdm;
    end
    save(['ppi_VOIGLM-' GLM_NAME '_results'  '.mat'], 'beta_mat', 'pval_mat', 'GLM_OPTS', 'voi_list');
end

%% plot results
if LOAD_GLM
    load(['ppi_VOIGLM-' GLM_NAME '_results'  '.mat'])
end

meanb_ML_mat = mean(beta_mat(:,1,:,:),1);
meanb_MH_mat = mean(beta_mat(:,2,:,:),1);
meanp_ML_mat = mean(pval_mat(:,1,:,:),1);
meanp_MH_mat = mean(pval_mat(:,2,:,:),1);
meanp_ML_mat(meanp_ML_mat > 0.05) = 1;
meanp_MH_mat(meanp_MH_mat > 0.05) = 1;

meanb_ML_mat = reshape(meanb_ML_mat, size(meanb_ML_mat,4), size(meanb_ML_mat,3));
meanb_MH_mat = reshape(meanb_MH_mat, size(meanb_MH_mat,4), size(meanb_MH_mat,3));
meanp_ML_mat = reshape(meanp_ML_mat, size(meanp_ML_mat,4), size(meanp_ML_mat,3));
meanp_MH_mat = reshape(meanp_MH_mat, size(meanp_MH_mat,4), size(meanp_MH_mat,3));

subplot(2,2,1)
imagesc(meanb_MH_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('Mean Beta - MSK HIGH');

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

subplot(2,2,4)
imagesc(meanp_ML_mat);
set(gca,'XTickLabel',voi_list)
set(gca,'YTickLabel',voi_list)
set(gca,'FontSize', 7)
set(gca,'XTick',1:length(voifile.VOI));
set(gca,'YTick',1:length(voifile.VOI));
set(gca,'XAxisLocation','top');
title('P Val - MSK LOW');








