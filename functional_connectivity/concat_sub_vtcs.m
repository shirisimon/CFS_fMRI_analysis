function concat_vtc = concat_sub_vtcs(subdir, INPUTPAT_VTC)
vtc = [];
vtcfiles = findFilesBVQX(subdir, INPUTPAT_VTC,struct('maxdepth',2));
for vt = 1:length(vtcfiles)
    vtc = BVQXfile(vtcfiles{vt});
    

end