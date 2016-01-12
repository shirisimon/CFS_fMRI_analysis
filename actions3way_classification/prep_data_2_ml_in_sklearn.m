
% import voi
clear all; close all; clc;

%% PARAMETERS
INPUTFILE_VOI = 'G:\study 3_CFS-fMRI_v2\data\mirror_loc_final.voi';
INPUTDIR      = 'G:\study 3_CFS-fMRI_v2\data';
INPUTPAT_VTC_TR = '*nmsk*mm.vtc'; 
INPUTPAT_VTC_TE = '*_msk*mm.vtc'; 
INPUTPAT_PRT_TR = '*nmsk*_acts*.prt'; 
INPUTPAT_PRT_TE = '*_msk*_acts*.prt'; 
VOIs2Clust    = 1:15; % excluding outliers based on PDIST results

voifile  = BVQXfile(INPUTFILE_VOI);
vtcfiles_TR = findFilesBVQX(INPUTDIR, INPUTPAT_VTC_TR, struct('maxdepth',3));
vtcfiles_TE = findFilesBVQX(INPUTDIR, INPUTPAT_VTC_TE, struct('maxdepth',3));
prtfiles_TR = findFilesBVQX(INPUTDIR, INPUTPAT_PRT_TR, struct('maxdepth',4));
prtfiles_TE = findFilesBVQX(INPUTDIR, INPUTPAT_PRT_TE, struct('maxdepth',4));
vtcfiles = [vtcfiles_TR; vtcfiles_TE];
prtfiles = [prtfiles_TR; prtfiles_TE];

%% LOAD VOI DATA:
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    for vo = 1:size(voifile.VOI,2)
        if vt==1; voilist{vo} = voifile.VOI(vo).Name; end
        voi_coords = tal2bv(voifile.VOI(vo).Voxels)';
        voi_mat = vtc.VOITimeCourseOrig(voi_coords, Inf);
        voi_mat = zscore(voi_mat{1});
        voidata{vt,vo} = voi_mat;
    end    
    clear voi_mat
    vtc.ClearObject();
end

%% LOAD PRT DATA
% conds2classify_msk  = {'msk_low_act1', 'msk_low_act2',  'msk_low_act3', ...
%                        'msk_high_act1','msk_high_act2','msk_high_act3'};
% conds2classify_nmsk = {'nmsk_low_act1', 'nmsk_low_act2',  'nmsk_low_act3', ...
%                        'nmsk_high_act1','nmsk_high_act2','nmsk_high_act3'};
conds2classify_msk  = {'msk_high_act1', 'msk_high_act2',  'msk_high_act3'};
conds2classify_nmsk = {'nmsk_high_act1', 'nmsk_high_act2',  'nmsk_high_act3'};
for p = 1:length(prtfiles)
    prt = BVQXfile(prtfiles{p});
    prt = prt.Cond;
    try conds_idx = getCondsTmps(prt, conds2classify_msk);
    catch
        conds_idx = getCondsTmps(prt, conds2classify_nmsk);
    end
    prtdata{p} = conds_idx; 
end

%%
% labels = [1:3, 1:3];
% train_test_idx = [1:44, 45:88];
save('data2skl_high_smoothed.mat', 'voidata', 'prtdata', 'vtcfiles', 'voilist', '-v7.3'); 

