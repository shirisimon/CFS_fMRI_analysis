function fileNames = findFiles(params)

rootDir = params.rootDir;
% define input files patterns and exclusion files patterns :
subjectPattern  = params.subjectPattern;
vtcFilePattern  = params.vtcFilePattern;
vtcFileExcPat   = params.vtcFileExcPat;
prtFilePattern  = params.prtFilePattern;
prtFileExcPat   = params.prtFileExcPat;
mskFilePattern  = params.mskFilePattern;
% define input folders and search depth:
subjectFolders  = findfiles(rootDir,subjectPattern,'dirs=1','depth=1');
vtcDepth = params.vtcDepth;
prtDepth = params.prtDepth;

%% find the files 
for i = 1:length(subjectFolders) % itirate over subjects to find files
    [~,subName]  = fileparts(subjectFolders{i});
    fileNames.(['s' subName]).folder = subjectFolders{i};
    rootDir = fileNames.(['s' subName]).folder;
    vtcFiles = findfiles(rootDir,vtcFilePattern,struct('maxdepth',vtcDepth) );
    prtFiles = findfiles(rootDir,prtFilePattern,struct('maxdepth',prtDepth));
    mskFiles = findfiles(rootDir,mskFilePattern);
    cnt = 1;
    for j = 1:length(vtcFiles)
        if isempty(regexp(vtcFiles{j},vtcFileExcPat)) % check if should exclude vtc
        fileNames.(['s' subName]).(['run' num2str(cnt)]).vtcFileName = vtcFiles{j};
        fileNames.(['s' subName]).(['run' num2str(cnt)]).prtFileName = prtFiles{j};
        fileNames.(['s' subName]).(['run' num2str(cnt)]).mskFileName = mskFiles{1}; % each run always has same mask
        cnt = cnt + 1;
        end
    end
end
%% 


%% print files found to a text file for inspection 
fid = fopen('filesFound.txt','w+');
subNames = fieldnames(fileNames);
fprintf(fid,'results from file walking: \n\n');
fprintf(fid,'found %d subjects:\n', length(subNames));

%print subject folders found at top of files
for i = 1:length(subNames)
    fprintf(fid,'%s\n',subNames{i});
end
fprintf(fid,'\n');

%print vtc prt pairs found at bottom of file
for j= 1:length(subNames)
    runNames = fieldnames(fileNames.(subNames{j}));
    for z = 1:length(runNames)
        if strcmp(runNames{z}(1:3),'run' ) % check if this is a run folder
            fprintf(fid,'subject: %s,  %s, vtc: \t\t%s\n',...
                subNames{j},runNames{z},fileNames.(subNames{j}).(runNames{z}).vtcFileName);
            fprintf(fid,'subject: %s,  %s, prt: \t\t%s\n',...
                subNames{j},runNames{z},fileNames.(subNames{j}).(runNames{z}).prtFileName);
            fprintf(fid,'subject: %s,  %s, msk: \t\t%s\n',...
                subNames{j},runNames{z},fileNames.(subNames{j}).(runNames{z}).mskFileName);
        end
        fprintf(fid,'\n');
    end
%     fprintf(fid,'\n');
end
fclose(fid); 
%% 

end