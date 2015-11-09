clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'D:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR      = 'D:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC  = '*_nmsk*mm.vtc'; % only msk vtcs
VOIs2Clust    = 1:15; % excluding outliers based on PDIST results

voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles = findFilesBVQX(INPUTDIR, INPUTPAT_VTC,struct('maxdepth',3));

% LOAD DATA:
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    for vo = 1:size(voifile.VOI,2)
        if vt==1
            voi_list{vo} = voifile.VOI(vo).Name;
        end
        voi_coords = tal2bv(voifile.VOI(vo).Voxels)';
        voi_rtc(:,vo) = zscore(vtc.VOITimeCourseOrig(voi_coords));
    end
    
    % for each subject:
    rdm(vt,:,:) = squareform(pdist(voi_rtc','correlation')); % corrcoef distance matrix
    
    clear voi_rtc
    vtc.ClearObject();
end
mrdm = mean(rdm);
mrdm = reshape(mrdm, size(mrdm,2), size(mrdm,3));


%% 1. PLOT PDIST (look at the distances):
[Y,e] = cmdscale(mrdm);
plot(Y(:,1),Y(:,2),'.');
text(Y(:,1),Y(:,2),voi_list);