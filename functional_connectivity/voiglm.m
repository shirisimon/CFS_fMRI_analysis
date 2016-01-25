function [beta, pval]  = voiglm(subdir, INPUTPAT_VTC, OUTPUTDIR_DMS, seed_name, contrasts, voi)

str_place = strsplit(subdir, '\'); % string place holder
SUB_NUM   = ['*' str_place{end}];
str_place = strfind(INPUTPAT_VTC, '*');
SDM_TYPE  = INPUTPAT_VTC(str_place(1):str_place(2)); 

vtcfiles = findFilesBVQX(subdir, INPUTPAT_VTC,struct('maxdepth',2));
sdmfiles = findFilesBVQX(OUTPUTDIR_DMS, [seed_name '_' SUB_NUM '*' SDM_TYPE]);

design_mat = [];
voi2_rtc   = [];
for vo2 = 1:size(voi.VOI,2);
    for vt = 1:length(vtcfiles)
        if vo2 == 1
            sdm = BVQXfile(sdmfiles{vt});
            design_mat = [design_mat; sdm.SDMMatrix];
        end
        vtc = BVQXfile(vtcfiles{vt});
        vtc.VTCData = zscore(vtc.VTCData);
        voi2_coords = tal2bv(voi.VOI(vo2).Voxels)';
        voi2_rtc    = [voi2_rtc; vtc.VOITimeCourseOrig(voi2_coords)];
        vtc.ClearObject; clear vtc;
        sdm.ClearObject; clear sdm;
    end
    X = design_mat;
    y = voi2_rtc;
    b = (X'*X)^(-1)*X'*y;
    e = y - X*b; 
    for c = 1:size(contrasts,2)
        t = (contrasts(:,c)'*b) / (sqrt(var(e)*contrasts(:,c)'*(X'*X)^(-1)*contrasts(:,c)));
        p = 2*(1-tcdf(abs(t),size(X,1) - size(X,2)));
        pval(vo2,c) = p; 
        beta(vo2,c) = b(contrasts(:,c)==1);
    end
    
end