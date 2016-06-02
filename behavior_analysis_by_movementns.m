function behavior_analysis_by_movementns()

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
        idx2remove = isnan(data(:,7));
        res_data = data;
        res_data(idx2remove,:) = [];
        res_data = response_mapping(res_data, map);
        res_data = assign_accuracy(res_data);
        
        low_opac = res_data(res_data(:,3) ~= 255,:);
        low_opac_act1 = low_opac(low_opac(:,1) == 0, :);
        low_opac_act2 = low_opac(low_opac(:,1) == 1, :);
        low_opac_act3 = low_opac(low_opac(:,1) == 2, :);
        % actions recognition:
        act1_rec_results(r,:) = rec_results(low_opac_act1);
        act2_rec_results(r,:) = rec_results(low_opac_act2);
        act3_rec_results(r,:) = rec_results(low_opac_act3);
        % actions perception confidence:
        act1_conf_results(r,:) = conf_results(low_opac_act1);
        act2_conf_results(r,:) = conf_results(low_opac_act2);
        act3_conf_results(r,:) = conf_results(low_opac_act3);
        
    end
    allsub_act1_rec(s,:) = mean(act1_rec_results);
    allsub_act2_rec(s,:) = mean(act2_rec_results);
    allsub_act3_rec(s,:) = mean(act3_rec_results);
    allsub_act1_conf(s,:) = mean(act1_conf_results);
    allsub_act2_conf(s,:) = mean(act2_conf_results);
    allsub_act3_conf(s,:) = mean(act3_conf_results);
    
end
end


function data = response_mapping(data, map)
for i = 1:size(data,1)
    switch data(i,7)
        case 2
            data(i,7) = map(2, data(i,9)+1);
        case 3
            data(i,7) = map(3, data(i,9)+1);
        case 4
            data(i,7) = map(4, data(i,9)+1);
        case 5
            data(i,7) = map(5, data(i,9)+1);
        case 6
            data(i,7) = map(6, data(i,9)+1);
    end
end
end


function data = assign_accuracy(data)
for i = 1:size(data,1)
    if data(i,7) == data(i,1)
        data(i,9) = 1;
    elseif data(i,7) == 3
        data(i,9) = 2;
    else
        data(i,9) = 0;
    end
end
end


function results = rec_results(data)
results = [sum(data(:,9)==1)/size(data,1), ...  % correct ans
    sum(data(:,9)==0)/size(data,1),  ... % wrong ans
    sum(data(:,9)==2)/size(data,1)]; ... % didn't see
end


function results = conf_results(data)
results = [sum(data(:,11)==0)/size(data,1), ...
    sum(data(:,11)==1)/size(data,1), ...
    sum(data(:,11)==2)/size(data,1), ...
    sum(data(:,11)==3)/size(data,1)];
end