%% addFoldsToSLinServer
clear all; close all; clc; 

files = dir('*.mat');
for s = 1:length(files)
    load(files(s).name);
    if ~exist('orglabels')
        orglabels = labels;
        labels(labels>3) = labels(labels>3) - 3;
        factor = [ones(1,size(labels,1)/2,1), 2*ones(1,size(labels,1)/2)];
        save(files(s).name, 'data', 'factor', 'labels', 'locations', 'map', ...
            'params', 'sName', 'vtcRes', '-v7.3');
    end
    clear orglabels
end