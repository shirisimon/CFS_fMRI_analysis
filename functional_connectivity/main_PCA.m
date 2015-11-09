clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR      = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_nmsk*mm.vtc'; % only msk vtcs
VOIs2Clust    = 1:11; % excluding outliers based on PDIST results

voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC,struct('maxdepth',3));

% LOAD DATA:
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    for vo = VOIs2Clust
        if vt==1
            voi_list{vo} = voifile.VOI(vo).Name;
        end
        voi_coords = tal2bv(voifile.VOI(vo).Voxels)';
        voi_rtc(:,vo) = zscore(vtc.VOITimeCourseOrig(voi_coords));
    end
    
    % for each subject:
    coefs{:,vt} = pca(voi_rtc(:,VOIs2Clust));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear voi_rtc
    vtc.ClearObject();
end
mkm = mean(km,2);
mkm(mkm>1.5) = 2;
mkm(mkm<1.5) = 1;
mkm_cell = [voi_list(VOIs2Clust)', num2cell(mkm)];


