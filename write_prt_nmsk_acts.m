function write_prt_nmsk_acts()
clear; close all; clc
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
subjects_dir_pattern = '389*';
files_pattern = '*_cnt2*.csv';
subjects_dir_depth = '1';
files_dir_depth = '2';

% params 2 prt
output_file_ext = 'nmsk2';
prt_type = 'vol';
convert2Vol = 1;
tr_length = 2;
tr_num = 338;
skip_tr = 2;
lag = 0;
predictors_names = {'nmsk_low_act1', 'nmsk_high_act1', ...
    'nmsk_low_act2', 'nmsk_high_act2', ...
    'nmsk_low_act3', 'nmsk_high_act3', ...
    'report'};
predictors_idx = [1 2 3]; % the importent colomns in logfile to prt


sub_folders = findfiles(rootdir, subjects_dir_pattern, 'dirs=1', ['depth=' subjects_dir_depth]);
for s = 1:length(sub_folders)
    datasets = findfiles(sub_folders{s}, files_pattern, ['depth=' files_dir_depth]);
    for r = 1:length(datasets)
        t = readtable(datasets{r});
        logfile_data = table2array(t);
        % fix csv data
        data = unique(logfile_data, 'rows', 'stable');
        idx2remove = isnan(data(:,1));
        data(idx2remove,:) = [];
        data = transform_tmps(data, skip_tr, lag, tr_length, prt_type);
        
        % define predictors:
        predictors = assign_tmps(data, predictors_names, predictors_idx);
        predictors(end+1).name = 'blnk';
        predictors(end).data = blnk_fill(predictors, tr_num,skip_tr);
        cons{1} = predictors(end).name; % add blnk data
        varibls_cons{1} = predictors(end).data;
        for p = 2:length(predictors)
            cons{p} = predictors(p-1).name;
            varibls_cons{p} = predictors(p-1).data;
        end
        
        % output fullfile
        path = [sub_folders{s}, '\' output_file_ext '\prt'];
        % path = [sub_folders{s}, '\loc\prt'];
        if convert2Vol; type = 'vol'; else type = 'ms'; end
        sub = strfind(sub_folders{s}, subjects_dir_pattern(1:2));
        sub = sub_folders{s}(sub:end);
        name = [sub '_' output_file_ext '_' type '_acts.prt'];
        ext = 'prt';
        if~exist(path, 'dir'); mkdir(path); end
        fid = fopen(fullfile(path, name), 'wt');
        colors = [ 0 0 0 ; 255 85 0 ; 0 170 127 ; 255 0 255 ; 170 255 127 ; 0 85 255 ; 255 0 0 ; 0 170 0] ;
        
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
        data(:,idx) = (round(data(:, idx)) / tr)+1 - skiptr + lag;
end
end


function data = assign_tmps(varibls, name, idx)
target_predictors = unique(varibls(:,idx(1)));
opacity_perdictors = unique(varibls(:,idx(2)));
c=0;
for t = 1:length(target_predictors)
    for m = 1:length(opacity_perdictors)
        c=c+1;
        data(c).name = name{c};
        curr_data = varibls(varibls(:,1) == target_predictors(t) & ...
            varibls(:,2) == opacity_perdictors(m) & ...
            varibls(:,3) ~= 1, 4);
        data(c).data(:,1) = unique(curr_data, 'stable');
        data(c).data(:,2) = data(c).data(:,1);
    end
end
data(c+1).name = 'catch';
data(c+1).data(:,1) = varibls(varibls(:,3) == 1, 4);
data(c+1).data(:,2) = data(c+1).data(:,1)+1;
% only in runs prt
% data(4).name = 'rec_response';
% data(4).data = varibls(:, [6 8]);
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
blnk(end+1,:) = [alltmp((end),2)+1 (trs-skiptr)];
blnk = sort(blnk,1);
end