ThisCatch = '406214';
load ModelName.mat
load FileName.mat

ResultsSummary = [];
for i = 1:numel(ModelName)

    ThisModel_MARRMoT = FileName{i};
    ThisModel = ModelName{i};
    try
    load(['EqualWeights/out/RRcal_out_AMALGAM_2obj_' ThisCatch '_' ThisModel '.mat']);
    AMALGAM_results = ResultsTable;
    load(['EqualWeights/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
    equal_CMAES_results = ResultsTable;
    load(['PrioritiseObj3/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
    priObj3_CMAES_results = ResultsTable;
    
    pp = [1, 1]; % perfect point
    ed = nan(size(AMALGAM_results,1), 1); % to store distance between pp and each point in the pareto front
    
    for ii = 1:size(AMALGAM_results,1)
        point = [table2array(AMALGAM_results(ii, 3)), table2array(AMALGAM_results(ii, 2))]; % make it a vector
        ed(ii) = EuclideanDistance(pp, point);
    end
    
    [~,ind] = min(ed);
    ctpp = [table2array(AMALGAM_results(ind, 3)), ...
            table2array(AMALGAM_results(ind, 2))]; % closest to perfect point (OFdry, OFnondry)
    
    % extract cmaes data point
    equal_point = [table2array(equal_CMAES_results(1, 4)), table2array(equal_CMAES_results(1, 3))];
    prior_point = [table2array(priObj3_CMAES_results(1, 4)), table2array(priObj3_CMAES_results(1, 3))];
    
    equal_of3 = table2array(equal_CMAES_results(:, 6));
    prior_of3 = table2array(priObj3_CMAES_results(:, 6));
    
    % distance between ctpp and of3 cal outcomes
    ed_equal = EuclideanDistance(ctpp, equal_point);
    ed_prior = EuclideanDistance(ctpp, prior_point);
    
    mean_ctpp = mean(ctpp);

    results = [real(equal_of3), real(prior_of3), real(ed_equal), real(ed_prior), real(mean(equal_point)), real(mean(prior_point))];
    ResultsSummary = [ResultsSummary; results];

    catch
        % If any error occurs, skip to the next iteration
        disp(['Model ' ThisModel_MARRMoT ' not available. Skipping...']);
        ResultsSummary = [ResultsSummary; [nan, nan, nan, nan, nan, nan]];
        continue;
    end
end

TradeoffsTable = array2table(ResultsSummary, "RowNames", ModelName, ...
    "VariableNames", {'equal_of3', 'prior_of3', 'ed_equal', 'ed_prior', 'equal_Q', 'prior_Q'});



% 
% % trade-offs plot
% figure;
% scatter(min(sqrt(2), ed_equal), mean_ctpp, [], equal_of3, 'filled', 'HandleVisibility', 'off')
% hold on;
% scatter(min(sqrt(2), ed_prior), mean_ctpp, [], prior_of3, 'filled', '^', 'HandleVisibility', 'off')
% xlim([0, sqrt(2)])
% ylim([0, 1])
% clim([0, 1])
% grid('on')
% c = colorbar;
% xlabel("Distance between CTPP and OF3 calibrated by CMA-ES")
% ylabel("Raw Q performance (i.e. mean of CTPP)")
% ylabel(c, "OF3")
% 
% % Create dummy scatter plots for the legend
% scatter(nan, nan, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'DisplayName', 'OF3 equal weight')
% scatter(nan, nan, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'DisplayName', 'OF3 higher weight')
% 
% % Add legend with only the marker shapes
% legend('show')



