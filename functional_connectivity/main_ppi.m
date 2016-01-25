
clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'G:\study3_CFS_fMRI_v2\data\mirror_loc_final_ppisorted.voi';
INPUTDIR      = 'G:\study3_CFS_fMRI_v2\data';
INPUTPAT_VTC  = '*_nmsk*mm.vtc'; % only msk vtcs
INPUTPAT_SDM  = '*_nmsk*mm*.sdm';
OUTPUTDIR_DMS = 'G:\study3_CFS_fMRI_v2\data\ppi_dms';
PRD_LIST      = '"seedvoi*NML" "seedvoi*NMH" "prdNML" "prdNMH" "seedvoi" "constant"';
GLM_NAME      = 'NM_manual_voisorted';
GENERATE_DMS  = 0; % generate sdms
GENERATE_GLM  = 1;
CONTRASTS     = [1 0 0 0 0 0; 0 1 0 0 0 0]; % contrasts for not-confound predictors

prdAcol       = 3; % colunm number in sdm file of condition A
prdBcol       = 4; % colunm number in sdm file of condition B
constant_col  = 7; % colunm number in sdm file of the constant

voifile  = BVQXfile(INPUTFILE_VOI);
subdirs  = findFilesBVQX(INPUTDIR, '3*',struct('dirs',1, 'maxdepth',1));
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC,struct('maxdepth',3));
sdmfiles = findFilesBVQX(INPUTDIR, INPUTPAT_SDM,struct('maxdepth',3));

%% MAIN
for vo1 = 1:size(voifile.VOI,2)
    %% 1. get seed voi rtc
    seed_coords = tal2bv(voifile.VOI(vo1).Voxels)';
    seed_name   = voifile.VOI(vo1).Name;
    fprintf('\n ppi processing for seed:  %s' , seed_name);
    
    if GENERATE_DMS
        for vt = 1:length(vtcfiles)
            vtc = BVQXfile(vtcfiles{vt});
            vtcnamepart = strsplit(vtcfiles{vt}, '\');
            vtcnamepart = strsplit(vtcnamepart{6}, '.');
            vtcnamepart = vtcnamepart{1};
            seed_rtc = zscore(vtc.VOITimeCourseOrig(seed_coords));
            
            %% 1. define interaction predictors
            sdm = BVQXfile(sdmfiles{vt});
            prdA = sdm.SDMMatrix(:,prdAcol);
            prdB = sdm.SDMMatrix(:,prdBcol);
            constant = sdm.SDMMatrix(:,constant_col);
            
            %% 2. define main predictors
            seedA = zscore(seed_rtc) .* prdA ;
            seedB = zscore(seed_rtc) .* prdB ;
            all_prds = [seedA, seedB, prdA, prdB, seed_rtc, constant];
            
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
    end
    %% 5. do manual VOIGLM and extract voi's betas
    if GENERATE_GLM
        for sub = 1:length(subdirs)
            [b, p]  = voiglm(subdirs{sub}, INPUTPAT_VTC, OUTPUTDIR_DMS, seed_name, CONTRASTS', voifile);
            beta_mat(sub,vo1,:,:) = b; % (vtc, voi1(seed), voi2, contrast)
            pval_mat(sub,vo1,:,:) = p;
        end
    end
    voi_list{vo1} = seed_name;
end
save(['ppi_VOIGLM-' GLM_NAME '.mat'], 'beta_mat', 'pval_mat', 'CONTRASTS', 'voi_list');






