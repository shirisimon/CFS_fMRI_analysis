function write_prt_nmsk()
clear all; close all; clc
% Input Paramaters:
% subject = int of subject number for output PRT filename
% run = string of run number for ouput PRT filename
% cons = cell array of strings with the names of regressors, should match
% (in terms of length):
% varibls = cell array of martices that has the start and end time of each
% event in msec in the format:
% [starttime1 , end time1 ;  starttime2 , endtime 2]
% colors = matrix of the colors you want each regressor to have in RGB
% CHANGE:
% file_pattern = loc\run
% tr_num
% predictors_names, predictors_col
% path and name to save
% more predictors in 'assign_tmps'


% params 2 find files
rootdir = 'D:\study 3_CFS-fMRI_v2\data\';
subjects_dir_pattern = '395*';
files_pattern = '*_cnt3*.csv';
subjects_dir_depth = '1';
files_dir_depth = '2';

% params 2 prt
output_file_ext = 'nmsk3';
predictors_idx = [30 0; 255 0];
prt_type = 'vol';
convert2Vol = 1;
tr_length = 2;
tr_num = 338;
skip_tr = 2;
lag = 0;
predictors_names = {'nmsk_low', 'nmsk_high', 'catch' };
predictors_cols = {'Opacity', 'Catch'};



sub_folders = findfiles(rootdir, subjects_dir_pattern, 'dirs=1', ['depth=' subjects_dir_depth]);
for s = 1:length(sub_folders)
    datasets = findfiles(sub_folders{s}, files_pattern, ['depth=' files_dir_depth]);
    for r = 1:length(datasets)
        logfile = importdata(datasets{r});
        % fix csv data
        data = unique(logfile.data, 'rows', 'stable');
        idx2remove = isnan(data(:,1));
        data(idx2remove,:) = [];
        data = transform_tmps(data, skip_tr, lag, tr_length, prt_type);
        %data(:,3) = data(:,3)+4000;
        
        % define predictors:
        % [~,predictors_idx] = ismember(predictors_cols, logfile.colheaders);
        predictors = assign_tmps(data, predictors_names, predictors_idx);
        predictors(end+1).name = 'blnk';
        predictors(end).data = blnk_fill(predictors, tr_num,skip_tr);
        cons{1} = predictors(end).name; % add blnk data
        cons{2} = 'msk_low'; predictors(1).name;
        cons{3} = 'msk_high';
        cons{4} = predictors(1).name;
        cons{5} = predictors(2).name;
        cons{6} = 'report';
        cons{7} = predictors(3).name;
        
        varibls_cons{1} = predictors(end).data;
        varibls_cons{2} = [];
        varibls_cons{3} = [];
        varibls_cons{4} = predictors(1).data;
        varibls_cons{5} = predictors(2).data;
        varibls_cons{6} = [];
        varibls_cons{7} = predictors(3).data;
        
        % output fullfile
        path = [sub_folders{s}, '\' output_file_ext '\prt'];
        % path = [sub_folders{s}, '\loc\prt'];
        if convert2Vol; type = 'vol'; else type = 'ms'; end
        sub = strfind(sub_folders{s}, subjects_dir_pattern(1:2));
        sub = sub_folders{s}(sub:end);
        name = [sub '_' output_file_ext '_' type '.prt'];
        ext = 'prt';
        if~exist(path, 'dir'); mkdir(path); end
        fid = fopen(fullfile(path, name), 'wt');
        colors = [ 0 	 0 	   0 ; 255 	 85 	 0 ; 0 	 170 	 127 ; 255 	 255 	 127 ; ...
                   85 	 170 	 255 ; 255 	 0 	 255 ; 255 	 170 	 127 ;  0   170   0] ;
        % write txt :
        fprintf(fid,'\n%s\n','FileVersion:    2') ;
        fprintf(fid,'\n%s\n','ResolutionOfTime:   vol     ');
        fprintf(fid,'\n%s\n','Experiment:         CFS_actions      ');
        fprintf(fid,'\n%s\n','BackgroundColor:    0 0 0     ');
        fprintf(fid,'\n%s\n','TextColor:          255 255 255    ');
        fprintf(fid,'\n%s\n','TimeCourseColor:    255 255 255      ');
        fprintf(fid,'\n%s\n','TimeCourseThick:    3    ');
        fprintf(fid,'\n%s\n','ReferenceFuncColor: 0 0 80      ');
        fprintf(fid,'\n%s\n','ReferenceFuncThick: 3      ');
        fprintf(fid,'\n%s\n',['NrOfConditions:   ' num2str(length(cons)) ]); % number of con / pred/ reg/ in PRT
        
        for i=1:length(cons)
            fprintf(fid,'\n%s\n', cons{i}) ;
            fprintf(fid,' %i  \n ', length(varibls_cons{i}));
            fprintf(fid,' %i \t %i \t  \n ',varibls_cons{i}');
            fprintf(fid,'%s %i \t %i \t %i  \n ','Color:	',colors(i,:));
        end
        
    end
end

end


function data = transform_tmps(data, skiptr, lag, tr, prt_type)
format short g;
idx = mod(sum(data,1),1);
idx = find(idx~=0);
switch prt_type
    case'ms'
        data(:,idx) = round((data(:, idx) - skiptr*tr) * ...
            1000) + lag*tr*1000  ;
    case'vol'
        data(:,idx) = (round(data(:, idx)) / tr)+1 - skiptr;
end
end


function data = assign_tmps(varibls, name, predictors)
%varibls(varibls(:,idx) == 3,idx) = 2;
% predictors = unique(varibls(:,idx));
% predictors = unique(varibls(:,idx),'rows');

for c = 1:size(predictors,1)
    data(c).name = name{c};
    data(c).data(:,1) = varibls(varibls(:,2) == predictors(c,1) & varibls(:,3) == predictors(c,2), 4); % for loc prt
    data(c).data(:,2) = data(c).data(:,1);
    % data(c).data = varibls(varibls(:,idx) == predictors(c), [3 4]); % for run prt
end
data(3).name = 'catch';
data(3).data(:,1) = varibls(varibls(:,3) == 1, 4);
data(3).data(:,2) = data(3).data(:,1)+2;
% data(5).name = 'conf_response';
% data(5).data = varibls(:, [8 10]);
end


function blnk = blnk_fill(predictors, trs, skiptr)
% set experiment's blank time points
alltmp = [];
for p = 1:length(predictors)
    alltmp = [alltmp; predictors(p).data];
end
alltmp = sort(alltmp, 1);

ind = 2 : length(alltmp);
blnk(1,:) = [1 alltmp(1,1)-1];
blnk(ind,:) = [alltmp((ind-1),2)+1 alltmp((ind),1)-1];
blnk(blnk(:,1) > blnk(:,2),:) = [];
blnk(end+1,:) = [alltmp((end),2)+1 (trs-skiptr)*2000];
blnk = sort(blnk,1);
end