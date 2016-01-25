function [condindices] = getCondsTmps(prtData, ...
    conditions2classify)


condindices = getCondIndices(prtData, conditions2classify);
% 
% extractedData = [];
% labels = [];
% for j = 1:length(condindices)
%     indices = condindices{j} + tmpFromOnset;
%     indices = indices(indices < size(vtcNorm,1)); % this is to make sure that the onsets in addition to the timepoint to exctract exist (like in a run that was cut)
%     indicesBaseline = condindices{j} + params.baselineToExtract;
%     indicesBaseline = indicesBaseline(indicesBaseline < size(vtcNorm,1)); % this is to make sure that the onsets in addition to the timepoint to exctract exist (like in a run that was cut)
%     % somtime basline and regular indices willl not match check that they
%     % match
%     [indices,indicesBaseline] = checkIndxVsBaseline(indices,indicesBaseline);
%     if params.extractPCT
%         dataTemp = vtcNorm(indices,:,:,:)./vtcNorm(indicesBaseline,:,:,:);
%     else
%         dataTemp = vtcNorm(indices,:,:,:);
%     end
%     extractedData = [extractedData  ; dataTemp];
%     labelsToAdd = ones(length(indices),1) * j;
%     labels = [labels ; labelsToAdd];
%     dataTemp = [] ; labelsToAdd = []; indices = []; indicesBaseline = [];
% end

end


function condindices = getCondIndices(cond_all, cond_sub)
for i = 1:length(cond_sub)
    cond = cond_sub{i};
    % find cond in cond_all
    for j = 1:length(cond_all)
        if strcmp(cond,cond_all(j).ConditionName)
            condindices{i} = cond_all(j).OnOffsets(:,1);
            break
        end
    end
end
end
