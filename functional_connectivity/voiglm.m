function [beta, pval]  = voiglm(sdm, vtc, contrasts, voi)

X   = sdm.SDMMatrix;

for vo2 = 1:size(voi.VOI,2);
    voi2_coords = tal2bv(voi.VOI(vo2).Voxels)';
    voi2_rtc    = zscore(vtc.VOITimeCourseOrig(voi2_coords));
    
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
