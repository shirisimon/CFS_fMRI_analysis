% params 2 find files
rootdir = 'G:\study3_CFS_fMRI_v2\data\';
subjects_dir_pattern = '3*';
files_pattern = '*_exp*.csv';
subjects_dir_depth = 1;
files_dir_depth = 2;


opts.dirs = 1; opts.depth = subjects_dir_depth;
sub_folders = findFilesBVQX(rootdir, subjects_dir_pattern, opts); % 'dirs=1', ['depth=' subjects_dir_depth]);
map = [0 1 2 3; 2 1 3 0; 1 3 0 2; 3 2 0 1; 0 3 2 1; 1 2 3 0];
for s = 1:length(sub_folders)
    opts.dirs = 0; opts.depth = files_dir_depth;
    datasets = findFilesBVQX(sub_folders{s}, files_pattern, opts);
    for r = 1:length(datasets)
        logfile = importdata(datasets{r});
        % fix csv data
        data = unique(logfile.data, 'rows', 'stable');
        idx2remove = isnan(data(:,1));
        data(idx2remove,:) = [];
        report_idx = data(:,4) == 0;
        data(report_idx,:) = [];
        no_report_trials(s,r) = sum(isnan(data(:,7)));
    end
end

resutls = sum(no_report_trials,2);