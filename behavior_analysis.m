function behavior_analysis()

% params 2 find files
rootdir = 'D:\study 3_CFS-fMRI_v2\data\';
subjects_dir_pattern = '3*';
files_pattern = '*_exp*.csv';
subjects_dir_depth = '1';
files_dir_depth = '2';

sub_folders = findfiles(rootdir, subjects_dir_pattern, 'dirs=1', ['depth=' subjects_dir_depth]);
map = [0 1 2 3; 2 1 3 0; 1 3 0 2; 3 2 0 1; 0 3 2 1; 1 2 3 0];
for s = 1:length(sub_folders)
    datasets = findfiles(sub_folders{s}, files_pattern, ['depth=' files_dir_depth]);
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
        
        high_opac = res_data(res_data(:,3) == 255,:);
        high_opac_rec_results(r,:) = [sum(high_opac(:,9)==1)/size(high_opac,1), ...  % correct ans
            sum(high_opac(:,9)==0)/size(high_opac,1), ... % wrong ans
            sum(high_opac(:,9)==2)/size(high_opac,1)];    % didnt see
        high_opac_conf_results(r,:) = [sum(high_opac(:,11)==0)/size(high_opac,1), ... % lowest confidence
            sum(high_opac(:,11)==1)/size(high_opac,1), ...
            sum(high_opac(:,11)==2)/size(high_opac,1), ...
            sum(high_opac(:,11)==3)/size(high_opac,1)]; % highest confidence 
        
        low_opac = res_data(res_data(:,3) ~= 255,:);
        low_opac_rec_results(r,:) = [sum(low_opac(:,9)==1)/size(low_opac,1), ...  % correct ans
            sum(low_opac(:,9)==0)/size(low_opac,1), ... % wrong ans
            sum(low_opac(:,9)==2)/size(low_opac,1)];    % didnt see
        low_opac_conf_results(r,:) = [sum(low_opac(:,11)==0)/size(low_opac,1), ...
            sum(low_opac(:,11)==1)/size(low_opac,1), ...
            sum(low_opac(:,11)==2)/size(low_opac,1), ...
            sum(low_opac(:,11)==3)/size(low_opac,1)];
    end
    allsub_rec_high_opac(s,:) = mean(high_opac_rec_results);
    allsub_rec_low_opac(s,:) = mean(low_opac_rec_results);
    allsub_conf_high_opac(s,:) = mean(high_opac_conf_results);
    allsub_conf_low_opac(s,:) = mean(low_opac_conf_results);
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