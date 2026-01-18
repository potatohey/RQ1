%% This is Ziqi's RQ1 final figures generation script
% including figures 5-12, and figures S1, S2, S4, S7, S8
% Date: 11/09/2025

%% add paths
addpath(genpath("Historic"))
addpath(genpath("Synthetic/02_StepChangeDesign"))
addpath(genpath("MARRMoT"))
exp_path = 'C:\Users\zizhang3\OneDrive\OneDrive - The University of Melbourne\Desktop\Ziqi_PhD\99_Documents\02_PhD papers\figs\';
%% save file switch
fileSave = 0; % default not save
%fileSave = 1;

%% Figure 5

% Load model info
load ModelName.mat
load FileName.mat

Catchments = {'406213', '406214', '407230'};
CatchColors = lines(3);
ModelName(45) = []; % delete prms

nModels = numel(ModelName);
nCatchments = numel(Catchments);

OF1_data = cell(nCatchments, 1);
OF2_data = cell(nCatchments, 1);
OF3_data = cell(nCatchments, 1);

for c = 1:nCatchments
    ThisCatch = Catchments{c};
    temp_OF1 = nan(nModels, 1);
    temp_OF2 = nan(nModels, 1);
    temp_OF3 = nan(nModels, 1);

    for i = 1:nModels
        ThisModel = ModelName{i};
        ThisModel_MARRMoT = FileName{i};

        load(['EqualWeights/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
        equal_CMAES_results = ResultsTable;

        temp_OF1(i) = real(table2array(equal_CMAES_results(1, 3)));
        temp_OF2(i) = real(table2array(equal_CMAES_results(1, 4)));
        temp_OF3(i) = real(table2array(equal_CMAES_results(1, 6)));

    end

    OF1_data{c} = temp_OF1;
    OF2_data{c} = temp_OF2;
    OF3_data{c} = temp_OF3;
end

% Combine into matrix
combined_data = [OF1_data; OF2_data; OF3_data]; % 9 cells
box_data = cell2mat(combined_data);
group_labels = repmat(1:9, nModels, 1); 
group_labels = group_labels(:);
box_data = box_data(:);
group_positions = [1 1.5 2.0, 3 3.5 4.0, 5 5.5 6.0]; % 3 groups, 3 catchments each
% Prepare colors per group (cycle per catchment)
group_colors = repmat(CatchColors, 3, 1); % 9 groups total

% Create boxplot with custom colors
figure;
h = boxplot(box_data, group_labels, ...
    'Positions', group_positions,...
    'Colors', group_colors, ...
    'Widths', 0.4); % <-- narrower boxes
set(h, {'LineWidth'}, {1.5});

% Apply color to each box individually
for j = 1:9
    patchHandles = findobj(gca, 'Tag', 'Box');
    % Set face color
    set(patchHandles(10-j), 'Color', group_colors(j, :));
end

% Remove default x-axis labels
set(gca, 'XTickLabel', []);

% Add custom labels centered under each OF group
xticks = get(gca, 'XTick');
mid_ticks = [mean(xticks(1:3)), mean(xticks(4:6)), mean(xticks(7:9))];
text(mid_ticks, repmat(-1.1,1,3), ...
    { {'OF1','Q (non-drought)'} , ...
      {'OF2','Q (drought)'} , ...
      {'OF3','multi-annual dynamics'} }, ...
    'HorizontalAlignment', 'center', 'FontSize', 11);

ylabel('Objective Function Value');
ylim([-1, 1.05])
title('OF1, OF2, OF3 across catchments', "FontSize", 12);
yticks([-1 -0.5 0 0.5 1]);
%grid on;
% Add vertical separator lines between OF groups
line([2.5 2.5], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off');
line([4.5 4.5], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off');
% Add legend for catchments
hold on;
for c = 1:nCatchments
    plot(NaN, NaN, 's', ...
        'MarkerFaceColor', "none", ...
        'MarkerEdgeColor', CatchColors(c,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', Catchments{c}, ...
        'MarkerSize', 8);
end
legend('Location', 'southeast', "FontSize", 10.5);

if fileSave
    exportgraphics(gcf, [exp_path 'Figure5.png']);
end

%% Figure 6

% Load Data
ThisCatch = '407230';
ThisModel_MARRMoT = 'm_17_penman_4p_3s';
ThisModel = 'penman';
load(['EqualWeights/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);

CMAES_theta = table2array(ResultsTable(1, 7:end));
timeseries = archive.CalData.timeseries;

% Set up MARRMoT
m = feval(ThisModel_MARRMoT);
m.input_climate = [timeseries.P,  timeseries.PET,  timeseries.PET]; % THIRD INPUT IS DUMMY VAR FOR NOW (SHOULD BE TEMP)
m.delta_t = 1;                
m.S0 = zeros([m.numStores, 1]);
m.theta = CMAES_theta; % Parameters
m.solver_opts; % initalise solver
my_solver_opts = m.default_solver_opts(); % set to default
m.solver_opts = my_solver_opts;

% run
[output_ex,...                                                             % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
 output_in,...                                                             % Internal model fluxes
 output_ss,...                                                             % Internal storages
 output_waterbalance] = ...                                                % Water balance check              
                        m.get_output();  

t = datetime(timeseries.year(:), timeseries.month(:), timeseries.day(:));

Storage = m.stores .* m.StoreSigns;
TotalStorage = sum(Storage, 2);
StoreLgd = ['TotalStorage', m.StoreNames{:}];

%
% Plot
figure;
set(gcf, 'Position', [100 100 1000 550], 'Color', 'w'); % Wider figure
ax1 = axes('Position', [0.02, 0.12, 0.30, 0.78]); 
imshow('C:\Users\zizhang3\OneDrive\OneDrive - The University of Melbourne\Desktop\Ziqi_PhD\99_Documents\02_PhD papers\figs\Penman.png');
title('Penman model structure', "FontSize", 15);
text(ax1, 0.02, -0.1, '(a)', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

ax2 = axes('Position', [0.40, 0.12, 0.58, 0.78]);

%figure;
h1 = plot(t, TotalStorage, "LineWidth", 2.8);
hold on;
h2 = plot(t, Storage(:, 1), "LineWidth", 1.4);
h3 = plot(t, Storage(:, 2), "LineWidth", 1.4);
h4 = plot(t, Storage(:, 3), "LineWidth", 1.4, "Color", 'black');

% Add gray vertical lines for pre-drought (95-97) and post-drought (08-10) periods
xline(datetime(1995, 1, 1), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xline(datetime(1997, 1, 1), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xline(datetime(2008, 1, 1), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xline(datetime(2010, 1, 1), '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1);
% Add text annotations
ylimits = ylim;
text(datetime(1996, 1, 1), ylimits(1) + (ylimits(2) - ylimits(1)) * 0.5, 'Pre-drought', 'Color', 'k', 'FontSize', 11, 'HorizontalAlignment', 'center');
text(datetime(2009, 1, 1), ylimits(1) + (ylimits(2) - ylimits(1)) * 0.5, 'Late drought', 'Color', 'k', 'FontSize', 11, 'HorizontalAlignment', 'center');
%Add arrows pointing to the text annotations
annotation('doublearrow', [0.435, 0.502], [0.51, 0.51], 'Color', 'k', 'LineWidth', 0.5);
annotation('doublearrow', [0.877, 0.945], [0.51, 0.51], 'Color', 'k', 'LineWidth', 0.5);

% Set x-axis limits
xlim([datetime(1994, 1, 1), datetime(2011, 1, 1)])

lgd = legend([h2 h3 h4 h1], {'S_r_z', 'S_d_e_f', 'C_r_e_s', 'Total Storage=S_r_z+S_d_e_f+C_r_e_s'}, "Location", "south");
lgd.FontSize = 11;
ax2.FontSize = 12;
ylabel('Storage in model (mm)', "FontSize", 13)
xlabel('Year', "FontSize", 13)
title([ThisCatch ' Penman'], "FontSize", 15)
text(ax2, 0.015, 0.057, '(b)', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
hold off;

if fileSave
    exportgraphics(gcf, [exp_path 'Figure6.png']);
end
%% Figure 7
% Catchments to plot
catchments = {'406213', '406214', '407230'};
labels = {'(a)', '(b)', '(c)'};
colors = get(groot,'defaultAxesColorOrder');

% Make a wide white-background figure
figure;
set(gcf, 'Position', [100 200 1200 380]);

% Create a tiled layout with compact spacing
t = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

for i = 1:length(catchments)
    ThisCatch = catchments{i};
    load(['Trade-offs/' ThisCatch, '_Tradeoffs.mat'])

    x = TradeoffsTable.equal_of3;
    y = TradeoffsTable.equal_Q;

    nexttile; hold on;
    scatter(x, y, 30, 'o', ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineWidth', 1.5);
  
    xlim([0, 1])
    ylim([-1, 1])
    
    xlabel('Multi-annual storage dynamics (OF3)')
    if i == 1
        ylabel('Q performance (Mean of OF1 & OF2)')
    end
    text(0.04, -0.7, labels{i}, 'FontSize', 16, 'FontWeight', 'bold');

    grid on
    set(gca, 'FontSize', 13)
    title([ThisCatch ' equal weights'], 'FontSize', 17)
end

if fileSave
    exportgraphics(gcf, [exp_path 'Figure7.png']);
end
%% Figure 8
ThisCatch = '406213';
load(['Trade-offs/' ThisCatch, '_Tradeoffs.mat'])

% Create figure
figure;
set(gcf, 'Position', [100 100 1400 550], 'Color', 'w'); % Wider figure

% LEFT plot
ax1 = subplot(1,2,1);
hold on;

% Data
x = TradeoffsTable.equal_of3;
y = TradeoffsTable.equal_Q;

x_thresh = 0.7;
y_thresh = 0.3;

% Group definitions
group3 = x < x_thresh & y > y_thresh;
group1 = x >= x_thresh & y > y_thresh;
group2 = y <= y_thresh;

% Shade top-left area
fill([0 x_thresh x_thresh 0], [y_thresh y_thresh 1 1], ...
    [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

% Scatter points
scatter(x(group1), y(group1), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x(group2), y(group2), 50, 's', 'filled', 'DisplayName', 'Fail');
scatter(x(group3), y(group3), 50, 'o', 'filled', 'DisplayName', 'Recalibration');

% Threshold lines
xline(x_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Labels
text(x_thresh + 0.02, 0.95, 'x = 0.7', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(0.03, y_thresh - 0.02, 'y = 0.3', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
title([ThisCatch ' equal weights (all models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 16)
grid on
box on
text(0.04, -0.9, '(a)', 'FontSize', 15.5, 'FontWeight', 'bold');

% RIGHT plot
ax2 = subplot(1,2,2);

% First create a grey background covering the whole right subplot
pos = get(ax2, 'Position');  % Get position of right subplot
expand = 0.055;               % Expand by 5% on all sides

annotation('rectangle', ...
    [pos(1)-expand pos(2)-1.9*expand pos(3)+1.5*expand pos(4)+3*expand], ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.15);

hold on;

% New data
x_prior = TradeoffsTable.prior_of3;
y_prior = TradeoffsTable.prior_Q;

x_thresh_prior = 0.7;
y_thresh_prior = 0.3;

x_g3 = x_prior(group3);
y_g3 = y_prior(group3);

% Define pass/fail
pass = x_g3 > x_thresh_prior & y_g3 > y_thresh_prior;
fail = ~pass;

% Scatter points
scatter(x_g3(pass), y_g3(pass), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x_g3(fail), y_g3(fail), 50, 's', 'filled', 'DisplayName', 'Fail');

% Threshold line
xline(x_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
%title([ThisCatch ' higher weights (32 models)'], 'FontSize', 18)
title(['Recalibration with higher weights (' num2str(length(x_g3)) ' models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 15.5)
grid on
box on
text(0.04, -0.9, '(b)', 'FontSize', 16, 'FontWeight', 'bold');

% Add a nice curved annotation arrow
annotation('textarrow', [0.36 0.515], [0.68 0.5], ...
    'String', '  Recalibration region', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'HeadLength', 12, 'HeadWidth', 12, ...
    'LineWidth', 2, 'Color', 'k');

if fileSave
    exportgraphics(gcf, [exp_path 'Figure8.png']);
end

%% Figure 9
load("RecoveryYears_subsPrior.mat")
RecoveryYearsTable_subsPrior(45,:) = []; % detlete prms
T = RecoveryYearsTable_subsPrior;   % 46x3 table of recovery years

catchmentLabels = string(T.Properties.VariableNames); 
if numel(catchmentLabels) ~= 3
    catchmentLabels = ["406213","406214","407230"];
end

A = table2array(T);

% --- BINNING --------------------------------------------------------------
edges = [0.5 1.5 2.5 3.5 inf];
binLabels = ["1 yr","2 yrs","3 yrs","\geq 4 yrs"];
countBins = @(x) histcounts(x, edges);

counts = zeros(3,4);
for c = 1:3
    counts(c,:) = countBins(A(:, c));
end

% --- PLOT -----------------------------------------------------------------
figure('Color','w', 'Position',[500 250 500 400]); hold on
xCenters  = 1:3;
barWidth  = 0.55;

% Define a colormap (4 shades of one color, e.g. blue)
shades = [
    0.85 0.92 1.00;   % 1 yr → very light blue
    0.55 0.75 0.98;   % 2 yrs → medium-light blue
    0.25 0.50 0.85;   % 3 yrs → mid-dark blue
    0.05 0.25 0.65];  % ≥4 yrs → darkest blue

b = bar(xCenters, counts, barWidth, 'stacked');

% Apply consistent colors and edges
for k = 1:4
    b(k).FaceColor = shades(k,:);
    b(k).EdgeColor = 'k';
    b(k).LineWidth = 1.2;
end

legend(b, binLabels, 'Location','eastoutside', ...
    'Orientation','vertical','Box','off','FontSize',11);

% Axes/labels
set(gca, 'XTick', xCenters, 'XTickLabel', catchmentLabels, 'FontSize', 12)
xlabel('Catchments','FontSize', 13)
ylabel('Number of models', 'FontSize', 13)
title('Years to 95% recovery', 'FontSize', 14)
xlim([0.5, 3.5])

% Y-limit sized to totals
totals = sum(counts,2);
ylim([0, max(totals)*1.15])

grid on; box on

if fileSave
    exportgraphics(gcf, [exp_path 'Figure9.png']);
end

%% Figure 10
% ========================================================================
%  LOAD DATA FOR PENMAN (Q_StepChange)
% ========================================================================
load("Q_StepChange.mat");
ThisCatch = 'x407230';   % used in Q_StepChange
modelName = 'penman';
paramSetting = 'equal';

modelData = Q_StepChange.(ThisCatch).(modelName);

% Extract Q_annual for Penman
Q_annual = modelData.(['Q_annual_' paramSetting]);
normalised_Q = (Q_annual - Q_annual(20)) / (Q_annual(30) - Q_annual(20));  

% ========================================================================
%  PENMAN MARRMoT SIMULATION
% ========================================================================
paramSetting = 'EqualWeights';
ThisModel = 'penman';
ThisModel_MARRMoT = 'm_17_penman_4p_3s';   % <- check your MARRMoT code list
ThisCatch_marrmot = '407230';

load(['StepChange' ThisCatch_marrmot '.mat']);
load([paramSetting '/out/RRcal_out_CMAES_' ThisCatch_marrmot '_' ThisModel '.mat']);
theta = table2array(ResultsTable(1, 7:end));

m = feval(ThisModel_MARRMoT);
m.input_climate = [StepChangeData.P.P, StepChangeData.PET.PET, StepChangeData.PET.PET]; % PET, dummy TEMP
m.delta_t = 1;                
m.S0 = zeros([m.numStores, 1]);
m.solver_opts = m.default_solver_opts();
m.theta = theta; 
    
% Run Penman
[output_ex, output_in, output_ss, output_waterbalance] = m.get_output();

Storage = m.stores .* m.StoreSigns;
TotalStorage = sum(Storage, 2);
StoreNames = [m.StoreNames{:}];

S = [StepChangeData.P(:, 1:3) array2table(Storage)];
S.Properties.VariableNames(4:end) = StoreNames;
S = [S array2table(TotalStorage)];
[~, lastRowIdx] = unique(S.Year, 'last');
S_annual = S(lastRowIdx, :);
S_annual.Month = [];
S_annual.Day = [];

% ========================================================================
%  GR4J RECOVERY FROM Q_StepChange
% ========================================================================
paramSetting_Q = 'equal';   % <- use lower-case 'equal' for Q_StepChange
modelName_gr4j = 'gr4j';
modelData_gr4j = Q_StepChange.(ThisCatch).(modelName_gr4j);

Q_annual_gr4j = modelData_gr4j.(['Q_annual_' paramSetting_Q]);
normalised_Q_gr4j = (Q_annual_gr4j - Q_annual_gr4j(20)) ./ ...
                    (Q_annual_gr4j(30) - Q_annual_gr4j(20));

% ========================================================================
%  GR4J MARRMoT SIMULATION
% ========================================================================
paramSetting = 'EqualWeights';  % <- only for loading calibration results
ThisModel_gr4j = 'gr4j';
ThisModel_MARRMoT_gr4j = 'm_07_gr4j_4p_2s';  % <- confirm correct identifier

load(['StepChange' ThisCatch_marrmot '.mat']);
load([paramSetting '/out/RRcal_out_CMAES_' ThisCatch_marrmot '_' ThisModel_gr4j '.mat']);
theta_gr4j = table2array(ResultsTable(1, 7:end));

mg = feval(ThisModel_MARRMoT_gr4j);
mg.input_climate = [StepChangeData.P.P, StepChangeData.PET.PET, StepChangeData.PET.PET]; 
mg.delta_t = 1;
mg.S0 = zeros([mg.numStores, 1]);
mg.solver_opts = mg.default_solver_opts();
mg.theta = theta_gr4j;

[output_ex_g, output_in_g, output_ss_g, output_waterbalance_g] = mg.get_output();

Storage_g = mg.stores .* mg.StoreSigns;
TotalStorage_g = sum(Storage_g, 2);
StoreNames_g = [mg.StoreNames{:}];

Sg = [StepChangeData.P(:, 1:3) array2table(Storage_g)];
Sg.Properties.VariableNames(4:end) = StoreNames_g;
Sg = [Sg array2table(TotalStorage_g)];
[~, lastRowIdx_g] = unique(Sg.Year, 'last');
S_annual_g = Sg(lastRowIdx_g, :);
S_annual_g.Month = [];
S_annual_g.Day = [];

% ========================================================================
%  FIGURE: GR4J (row 1) and PENMAN (row 2)
% ========================================================================
f = figure;
f.Position = [300 150 1150 700];   % wider + taller

% --- Panel positions: [left bottom width height]
pos_a = [0.07 0.58 0.2 0.35];   % (a) GR4J Recovery
pos_b = [0.35 0.58 0.62 0.35];   % (b) GR4J Storage
pos_c = [0.07 0.1 0.2 0.35];   % (c) Penman Recovery
pos_d = [0.35 0.1 0.62 0.35];   % (d) Penman Storage

tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% ---------- (a) GR4J Recovery ----------
ax1 = axes('Position',pos_a);
plot(20:30, normalised_Q_gr4j(20:30)*100, '-o', 'Color', [0.47 0.67 0.19], ...
    'LineWidth', 2, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerSize', 4);
grid on; xlim([20 30]); ylim([0 100]);
ax1.FontSize = 12; ax1.XTick = 20:2:30; ax1.XTickLabel = string(20:2:30);
xlabel('Synthetic Year', "FontSize", 14);
ylabel('Recovery of Q (%)', "FontSize", 14);
title('Streamflow recovery (GR4J)', 'FontSize', 14);
text(ax1, 0.02, 0.1, '(a)', 'Units','normalized','FontSize',14,'FontWeight','bold');

% ---------- (b) GR4J Storage ----------
ax2 = axes('Position',pos_b);
hg4 = plot(TotalStorage_g, "LineWidth", 2, "DisplayName", "Total Storage"); hold on;
hg1 = plot(Storage_g(:,1), "LineWidth", 1.1, "DisplayName", StoreNames_g(1));
hg2 = plot(Storage_g(:,2), "LineWidth", 1.1, "DisplayName", StoreNames_g(2));
grid on;
xtick_data   = linspace(0,14600,9); xtick_labels = 0:5:40;
set(ax2,'XTick',xtick_data,'XTickLabel',xtick_labels)
ax2.FontSize = 12;
ylabel('Storage in Model (mm)', "FontSize", 14);
xlabel('Synthetic Year', "FontSize", 14);
title(sprintf('Simulated storage (GR4J %s)', ThisCatch_marrmot), "FontSize", 14);
text(ax2, 0.02, 0.1, '(b)', 'Units','normalized','FontSize',14,'FontWeight','bold');

% Legend GR4J

legend([hg1 hg2 hg4], {'S store', 'R store', 'Total Storage = S + R'}, ...
        "Location","northeast","FontSize",11);


% ---------- (c) Penman Recovery ----------
ax3 = axes('Position',pos_c);
plot(20:30, normalised_Q(20:30)*100, '-o', 'Color', [0.47 0.67 0.19], ...
    'LineWidth', 2, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerSize', 4);
grid on; xlim([20 30]); ylim([0 100]);
ax3.FontSize = 12; ax3.XTick = 20:2:30; ax3.XTickLabel = string(20:2:30);
xlabel('Synthetic Year', "FontSize", 14);
ylabel('Recovery of Q (%)', "FontSize", 14);
title('Streamflow recovery (Penman)', 'FontSize', 14);
text(ax3, 0.02, 0.1, '(c)', 'Units','normalized','FontSize',14,'FontWeight','bold');

% ---------- (d) Penman Storage ----------
ax4 = axes('Position',pos_d);
h4 = plot(TotalStorage, "LineWidth", 2, "DisplayName", "Total Storage"); hold on;
h1 = plot(Storage(:,1), "LineWidth", 1.1, "DisplayName", StoreNames(1));
h2 = plot(Storage(:,2), "LineWidth", 1.1, "DisplayName", StoreNames(2));
h3 = plot(Storage(:,3), "LineWidth", 1.1, "Color","k", "DisplayName", StoreNames(3));
legend([h1 h2 h3 h4], {'S_r_z', 'S_d_e_f', 'C_r_e_s', 'Total Storage=S_r_z+S_d_e_f+C_r_e_s'}, "Location", "southeast", "FontSize", 11);

grid on;
set(ax4,'XTick',xtick_data,'XTickLabel',xtick_labels)
ax4.FontSize = 12;
ylabel('Storage in Model (mm)', "FontSize", 14);
xlabel('Synthetic Year', "FontSize", 14);
title(sprintf('Simulated storage (Penman %s)', ThisCatch_marrmot), "FontSize", 14);
text(ax4, 0.02, 0.1, '(d)', 'Units','normalized','FontSize',14,'FontWeight','bold');

if fileSave
    exportgraphics(gcf, [exp_path 'Figure10.png']);
end
%% Figure 11
% Model counts
with_count = 9;
without_count = 37;
start_angle = 90;
angle_with = 360 * with_count / (with_count + without_count);
angle_without = 360 - angle_with;

% Test data: [with_hist, with_syn, without_hist, without_syn]
test_data = [
    8, 4, 5, 1;
    6, 3, 1, 0;
    8, 5, 2, 1
];

% Colors
blue = [201 223 229]/255;
grey = [210 210 210]/255;

% Shifts
x_shift = 0.02;
y_shift = -0.02;

% Figure
f = figure;
f.Position = [600 300 1000 520];
annotation('line', [0.235 0.235] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
annotation('line', [0.475 0.475] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
annotation('line', [0.727 0.727] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);

annotation('textbox', [0.257, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '406213', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.505, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '406214', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.755, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '407230', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.015, 0.22, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '(HSC = hypothesized structual components)', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 12);

%-----------------------No. models----------------------
annotation('textbox', [0, 0.68, 0.2, 0.05], ...
    'String', '9', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', [0.4, 0.6, 1], 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.05, 0.5, 0.2, 0.05], ...
    'String', '37', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', [130 130 130]/255, 'FontWeight','bold',...
    'FontSize', 15);

% 406213
annotation('textbox', [0.24, 0.85, 0.2, 0.05], ...
    'String', num2str(test_data(1,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.333, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(1,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.176, 0.27-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(1,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.23, 0.44-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(1,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

% 406214
annotation('textbox', [0.5, 0.858, 0.2, 0.05], ...
    'String', num2str(test_data(2,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.582, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(2,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.412, 0.385-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(2,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.48, 0.44-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(2,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);
%407230
annotation('textbox', [0.727, 0.827, 0.2, 0.05], ...
    'String', num2str(test_data(3,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.833, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(3,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.73, 0.395, 0.2, 0.05], ...
    'String', num2str(test_data(3,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.665, 0.29, 0.2, 0.05], ...
    'String', num2str(test_data(3,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

% -------- Combined pie subplot (left-middle) --------
combined_ax = axes('Position', [0.015 + x_shift, 0.4 + y_shift, 0.2, 0.45]);
hold(combined_ax, 'on'); axis(combined_ax, 'equal'); axis(combined_ax, 'off');

offset = 1.5; % for visual separation
fill_arc(0, 0, 1, start_angle + angle_with, start_angle + 360, grey, 'none'); % grey base
fill_arc(-0.07, 0.07, 1, start_angle - offset, start_angle + angle_with - offset, blue, 'none'); % blue offset

p1 = patch(NaN, NaN, [201 223 229]/255, 'EdgeColor', 'none'); % blue
p2 = patch(NaN, NaN, [210 210 210]/255, 'EdgeColor', 'none'); % grey

% Add legend
legend([p1 p2], {'include all HSC', 'not include all HSC'}, ...
    'Box', 'off', ...
    'Location', [0.055 + x_shift+0.01, 0.3 + y_shift, 0.1, 0.1], ...
    'FontSize', 13, 'FontWeight','bold');

% -------- Main 3-column layout --------
full_width = 0.2;
gap = 0.05;
x_start = 0.35 + x_shift;

for i = 1:3
    d = test_data(i, :);
    wh = d(1); ws = d(2); wth = d(3); wts = d(4);

    x_center = x_start + (i-1)*(full_width + gap);
    left_with = x_center - 0.2/4;
    left_without = x_center - 0.2/2;

    % --- Upper subplot (blue) ---
    ax1 = axes('Position', [left_with, 0.7 + y_shift, 0.2/2, 0.45/2]);
    hold(ax1, 'on'); axis(ax1, 'equal'); axis(ax1, 'off');
    fill_arc(0, 0, 1, start_angle, start_angle + angle_with, blue, 'none');
    fill_arc(0, 0, 1.05, start_angle, start_angle + angle_with * (wh / with_count), 'none', 'k');
    fill_arc(0, 0, 1.1, start_angle, start_angle + angle_with * (ws / with_count), 'none', 'r');

    % --- Lower subplot (grey) ---
    ax2 = axes('Position', [left_without, 0.19 + y_shift-0.04, 0.2, 0.45]);
    hold(ax2, 'on'); axis(ax2, 'equal'); axis(ax2, 'off');
    fill_arc(0, 0, 1, start_angle + angle_with, 360 + start_angle, grey, 'none');
    fill_arc(0, 0, 1.05, start_angle + angle_with, start_angle + angle_with + angle_without * (wth / without_count), 'none', 'k');
    if i == 3
        fill_arc(0, 0, 1.1, 1.3+start_angle + angle_with + angle_without * (wth / without_count), ...
            start_angle + angle_with + angle_without * (wth / without_count) + angle_without * (wts / without_count), 'none', 'r');
    else
        fill_arc(0, 0, 1.1, start_angle + angle_with, start_angle + angle_with + angle_without * (wts / without_count), 'none', 'r');
    end
end

% Create square patches with no fill and only colored edges
p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'k', 'LineWidth', 1.2); % black edge only
p_spacer = patch(NaN, NaN, 'w', 'FaceColor', 'none', 'EdgeColor', 'none', 'Visible', 'off'); % invisible spacer
p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'r', 'LineWidth', 1.2); % red edge only

% Add legend with two rows
legend([p1 p_spacer p2], ...
       {'Passed historic test', '', 'Passed synthetic test'}, ...
       'Box', 'off', ...
       'FontSize', 14, ...
       'FontWeight','bold',...
       'NumColumns', 3, ...
       'Position', [0.35 + x_shift, 0.12 + y_shift-0.02, 0.5, 0.03]);

if fileSave
    exportgraphics(gcf, [exp_path 'Figure11.png']);
end

% % ---------- Helper Function ----------
% function fill_arc(xc, yc, r, theta1, theta2, facecolor, edgecolor)
%     t = linspace(deg2rad(theta1), deg2rad(theta2), 100);
%     x = [xc + r * cos(t), xc];
%     y = [yc + r * sin(t), yc];
%     if ischar(facecolor) && strcmp(facecolor, 'none')
%         fill(x, y, 'w', 'FaceColor', 'none', 'EdgeColor', edgecolor, 'LineWidth', 1.2);
%     else
%         fill(x, y, facecolor, 'EdgeColor', edgecolor, 'LineWidth', 1.2);
%     end
% end

%% Figure 12
% ----------- First Catchment: 406214 ----------- %%
load("Q_StepChange.mat");
ThisCatch = 'x406214';
modelName = 'lascam';
paramSetting = 'prior';

modelData = Q_StepChange.(ThisCatch).(modelName);
Q_annual = modelData.(['Q_annual_' paramSetting]);
normalised_Q1 = (Q_annual - Q_annual(20)) / (Q_annual(30) - Q_annual(20));

paramSetting = 'PrioritiseObj3';
ThisModel = 'lascam';
ThisModel_MARRMoT = 'm_23_lascam_24p_3s';
ThisCatch = '406214';
load(['StepChange' ThisCatch '.mat']);
load([paramSetting '/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
theta = table2array(ResultsTable(1, 7:end));

m1 = feval(ThisModel_MARRMoT);
m1.input_climate = [StepChangeData.P.P, StepChangeData.PET.PET, StepChangeData.PET.PET];
m1.delta_t = 1;
m1.S0 = zeros([m1.numStores, 1]);
m1.solver_opts = m1.default_solver_opts();
m1.theta = theta;

[~, ~, ~, ~] = m1.get_output();
Storage1 = m1.stores .* m1.StoreSigns;
TotalStorage1 = sum(Storage1, 2);

S = [StepChangeData.P(:, 1:3), array2table(Storage1)];
StoreNames1 = [m1.StoreNames{:}];
S.Properties.VariableNames(4:end) = StoreNames1;
S = [S, array2table(TotalStorage1)];
S.Properties.VariableNames{end} = 'TotalStorage';

[~, lastRowIdx] = unique(S.Year, 'last');
S_annual1 = S(lastRowIdx, :);
S_annual1.Month = [];
S_annual1.Day = [];

% ----------- Second Catchment: 407230 ----------- %%
load("Q_StepChange.mat");
ThisCatch = 'x407230';
modelName = 'lascam';
paramSetting = 'equal';

modelData = Q_StepChange.(ThisCatch).(modelName);
Q_annual = modelData.(['Q_annual_' paramSetting]);
normalised_Q2 = (Q_annual - Q_annual(20)) / (Q_annual(30) - Q_annual(20));

paramSetting = 'EqualWeights';
ThisModel = 'lascam';
ThisModel_MARRMoT = 'm_23_lascam_24p_3s';
ThisCatch = '407230';
load(['StepChange' ThisCatch '.mat']);
load([paramSetting '/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
theta = table2array(ResultsTable(1, 7:end));

m2 = feval(ThisModel_MARRMoT);
m2.input_climate = [StepChangeData.P.P, StepChangeData.PET.PET, StepChangeData.PET.PET];
m2.delta_t = 1;
m2.S0 = zeros([m2.numStores, 1]);
m2.solver_opts = m2.default_solver_opts();
m2.theta = theta;

[~, ~, ~, ~] = m2.get_output();
Storage2 = m2.stores .* m2.StoreSigns;
TotalStorage2 = sum(Storage2, 2);

S = [StepChangeData.P(:, 1:3), array2table(Storage2)];
StoreNames2 = [m2.StoreNames{:}];
S.Properties.VariableNames(4:end) = StoreNames2;
S = [S, array2table(TotalStorage2)];
S.Properties.VariableNames{end} = 'TotalStorage';

[~, lastRowIdx] = unique(S.Year, 'last');
S_annual2 = S(lastRowIdx, :);
S_annual2.Month = [];
S_annual2.Day = [];

% ----------- Create Custom Figure with Image and Two Storage Panels ----------- %%
f = figure;
f.Position = [600 200 900 750];
f.Color = 'w';

% Panel positions
left_panel_img = [-0.02 0.22 0.4 0.6];     % (a): spans two rows
top_right = [0.45 0.58 0.5 0.35];         % (b): top-right
bot_right = [0.45 0.08 0.5 0.35];         % (c): bottom-right

% --------- (a) Image Panel ---------
ax_img = axes(f, 'Position', left_panel_img);
imshow('C:\Users\zizhang3\OneDrive\OneDrive - The University of Melbourne\Desktop\Ziqi_PhD\99_Documents\02_PhD papers\figs\LASCAM.png');
axis off
text(ax_img, 0.02, -0.05, '(a)', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');
title('LASCAM model structure', 'FontSize', 15, 'FontWeight', 'bold')

ax2 = axes(f, 'Position', top_right);
h4 = plot(TotalStorage1, "LineWidth", 2, "DisplayName", "Total Storage");
hold on;
h1 = plot(Storage1(:, 1), "LineWidth", 1.1, "DisplayName", "F store");
h2 = plot(Storage1(:, 2), "LineWidth", 1.1, "DisplayName", "A store");
h3 = plot(Storage1(:, 3), "LineWidth", 1.1, "Color", "k", "DisplayName", "B store");
hold off;
grid on;
xtick_data   = linspace(0,14600,9); xtick_labels = 0:5:40;
set(ax2,'XTick',xtick_data,'XTickLabel',xtick_labels)
ax2.FontSize = 12;
xlabel('Synthetic Year', "FontSize", 14);
ylabel('Storage in Model (mm)', "FontSize", 14);
ylim([0, 1600]);
title('Simulated storage components (406214)', "FontSize", 15);
legend([h1 h2 h3 h4], {'F store', 'A store', 'B store', 'Total Storage = F+A+B'}, "Location", "east", "FontSize", 12);
text(ax2, -0.135, -0.13, '(b)', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

% --------- (c) Storage - 407230 ---------
ax3 = axes(f, 'Position', bot_right);
plot(TotalStorage2, "LineWidth", 2, "DisplayName", "Total Storage");
hold on;
plot(Storage2(:, 1), "LineWidth", 1.1, "DisplayName", "F store");
plot(Storage2(:, 2), "LineWidth", 1.1, "DisplayName", "A store");
plot(Storage2(:, 3), "LineWidth", 1.1, "Color", "k", "DisplayName", "B store");
hold off;
grid on;
xtick_data   = linspace(0,14600,9); xtick_labels = 0:5:40;
set(ax3,'XTick',xtick_data,'XTickLabel',xtick_labels)
ax3.FontSize = 12;
xlabel('Synthetic Year', "FontSize", 14);
ylabel('Storage in Model (mm)', "FontSize", 14);
title('Simulated storage components (407230)', "FontSize", 15);
text(ax3, -0.135, -0.13, '(c)', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold');

if fileSave
    exportgraphics(gcf, [exp_path 'Figure12.png']);
end
%% Figure S1
ThisCatch = '406214';
load(['Trade-offs/' ThisCatch, '_Tradeoffs.mat'])

% Create figure
figure;
set(gcf, 'Position', [100 100 1400 550], 'Color', 'w'); % Wider figure

% LEFT plot
ax1 = subplot(1,2,1);
hold on;

% Data
x = TradeoffsTable.equal_of3;
y = TradeoffsTable.equal_Q;

x_thresh = 0.7;
y_thresh = 0.3;

% Group definitions
group3 = x < x_thresh & y > y_thresh;
group1 = x >= x_thresh & y > y_thresh;
group2 = y <= y_thresh;

% Shade top-left area
fill([0 x_thresh x_thresh 0], [y_thresh y_thresh 1 1], ...
    [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

% Scatter points
scatter(x(group1), y(group1), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x(group2), y(group2), 50, 's', 'filled', 'DisplayName', 'Fail');
scatter(x(group3), y(group3), 50, 'o', 'filled', 'DisplayName', 'Recalibration');

% Threshold lines
xline(x_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Labels
text(x_thresh + 0.02, 0.95, 'x = 0.7', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(0.03, y_thresh - 0.02, 'y = 0.3', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
title([ThisCatch ' equal weights (all models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 16)
grid on
box on
text(0.04, -0.9, '(a)', 'FontSize', 15.5, 'FontWeight', 'bold');

% RIGHT plot
ax2 = subplot(1,2,2);

% First create a grey background covering the whole right subplot
pos = get(ax2, 'Position');  % Get position of right subplot
expand = 0.055;               % Expand by 5% on all sides

annotation('rectangle', ...
    [pos(1)-expand pos(2)-1.9*expand pos(3)+1.5*expand pos(4)+3*expand], ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.15);

hold on;

% New data
x_prior = TradeoffsTable.prior_of3;
y_prior = TradeoffsTable.prior_Q;

x_thresh_prior = 0.7;
y_thresh_prior = 0.3;

x_g3 = x_prior(group3);
y_g3 = y_prior(group3);

% Define pass/fail
pass = x_g3 > x_thresh_prior & y_g3 > y_thresh_prior;
fail = ~pass;

% Scatter points
scatter(x_g3(pass), y_g3(pass), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x_g3(fail), y_g3(fail), 50, 's', 'filled', 'DisplayName', 'Fail');

% Threshold line
xline(x_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
%title([ThisCatch ' higher weights (32 models)'], 'FontSize', 18)
title(['Recalibration with higher weights (' num2str(length(x_g3)) ' models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 15.5)
grid on
box on
text(0.04, -0.9, '(b)', 'FontSize', 16, 'FontWeight', 'bold');

% Add a nice curved annotation arrow
annotation('textarrow', [0.36 0.515], [0.68 0.5], ...
    'String', '  Recalibration region', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'HeadLength', 12, 'HeadWidth', 12, ...
    'LineWidth', 2, 'Color', 'k');

if fileSave
    exportgraphics(gcf, [exp_path 'FigureS1.png']);
end

%% Figure S2
ThisCatch = '407230';
load(['Trade-offs/' ThisCatch, '_Tradeoffs.mat'])

% Create figure
figure;
set(gcf, 'Position', [100 100 1400 550], 'Color', 'w'); % Wider figure

% LEFT plot
ax1 = subplot(1,2,1);
hold on;

% Data
x = TradeoffsTable.equal_of3;
y = TradeoffsTable.equal_Q;

x_thresh = 0.7;
y_thresh = 0.3;

% Group definitions
group3 = x < x_thresh & y > y_thresh;
group1 = x >= x_thresh & y > y_thresh;
group2 = y <= y_thresh;

% Shade top-left area
fill([0 x_thresh x_thresh 0], [y_thresh y_thresh 1 1], ...
    [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');

% Scatter points
scatter(x(group1), y(group1), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x(group2), y(group2), 50, 's', 'filled', 'DisplayName', 'Fail');
scatter(x(group3), y(group3), 50, 'o', 'filled', 'DisplayName', 'Recalibration');

% Threshold lines
xline(x_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Labels
text(x_thresh + 0.02, 0.95, 'x = 0.7', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(0.03, y_thresh - 0.02, 'y = 0.3', ...
    'FontSize', 14, 'FontWeight', 'normal', ...
    'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
title([ThisCatch ' equal weights (all models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 16)
grid on
box on
text(0.04, -0.9, '(a)', 'FontSize', 15.5, 'FontWeight', 'bold');

% RIGHT plot
ax2 = subplot(1,2,2);

% First create a grey background covering the whole right subplot
pos = get(ax2, 'Position');  % Get position of right subplot
expand = 0.055;               % Expand by 5% on all sides

annotation('rectangle', ...
    [pos(1)-expand pos(2)-1.9*expand pos(3)+1.5*expand pos(4)+3*expand], ...
    'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.15);

hold on;

% New data
x_prior = TradeoffsTable.prior_of3;
y_prior = TradeoffsTable.prior_Q;

x_thresh_prior = 0.7;
y_thresh_prior = 0.3;

x_g3 = x_prior(group3);
y_g3 = y_prior(group3);

% Define pass/fail
pass = x_g3 > x_thresh_prior & y_g3 > y_thresh_prior;
fail = ~pass;

% Scatter points
scatter(x_g3(pass), y_g3(pass), 50, '^', 'filled', 'DisplayName', 'Pass');
scatter(x_g3(fail), y_g3(fail), 50, 's', 'filled', 'DisplayName', 'Fail');

% Threshold line
xline(x_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
yline(y_thresh_prior, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

xlim([0, 1])
ylim([-1, 1])
set(gca, 'FontSize', 14)
xlabel('Multi-annual storage dynamics (OF3)', 'FontSize', 16)
ylabel('Q performance (Mean of OF1 & OF2)', 'FontSize', 16)
%title([ThisCatch ' higher weights (32 models)'], 'FontSize', 18)
title(['Recalibration with higher weights (' num2str(length(x_g3)) ' models)'], 'FontSize', 16.5)
yticks([-1 -0.5 0 0.5 1]);
legend('Location', 'southeast', 'FontSize', 15.5)
grid on
box on
text(0.04, -0.9, '(b)', 'FontSize', 16, 'FontWeight', 'bold');

% Add a nice curved annotation arrow
annotation('textarrow', [0.36 0.515], [0.68 0.5], ...
    'String', '  Recalibration region', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'HeadLength', 12, 'HeadWidth', 12, ...
    'LineWidth', 2, 'Color', 'k');

if fileSave
    exportgraphics(gcf, [exp_path 'FigureS2.png']);
end

%% Figure S4
ThisCatch = '406213';
ThisModel_MARRMoT = 'm_19_australia_8p_3s';
ThisModel = 'australia';
load(['PrioritiseObj3/out/RRcal_out_CMAES_' ThisCatch '_' ThisModel '.mat']);
%
CMAES_theta = table2array(ResultsTable(1, 7:end));
timeseries = archive.CalData.timeseries;

% Set up MARRMoT
m = feval(ThisModel_MARRMoT);
m.input_climate = [timeseries.P,  timeseries.PET,  timeseries.PET]; % THIRD INPUT IS DUMMY VAR FOR NOW (SHOULD BE TEMP)
m.delta_t = 1;                
m.S0 = zeros([m.numStores, 1]);
m.theta = CMAES_theta; % Parameters
m.solver_opts; % initalise solver
my_solver_opts = m.default_solver_opts(); % set to default
m.solver_opts = my_solver_opts;

% run
[output_ex,...                                                             % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
 output_in,...                                                             % Internal model fluxes
 output_ss,...                                                             % Internal storages
 output_waterbalance] = ...                                                % Water balance check              
                        m.get_output();  

t = datetime(timeseries.year(:), timeseries.month(:), timeseries.day(:));

Storage = m.stores .* m.StoreSigns;
TotalStorage = sum(Storage, 2);
StoreLgd = ['TotalStorage', m.StoreNames{:}];
%
yearMonth = dateshift(t, 'start', 'month');

% Aggregate Q values monthly
[uniqueMonths, ~, monthIdx] = unique(yearMonth);
monthlyQobs = accumarray(monthIdx, timeseries.Q, [], @sum);
monthlyQsim = accumarray(monthIdx, output_ex.Q, [], @sum);

%
% Create figure
figure;
set(gcf, 'Position', [100 100 1200 550], 'Color', 'w'); % Wider figure
ax1 = axes('Position', [-0.02, 0.08, 0.3, 0.82]); 
imshow('C:\Users\zizhang3\OneDrive\OneDrive - The University of Melbourne\Desktop\Ziqi_PhD\99_Documents\02_PhD papers\figs\Australia.png');
title('Australia model structure', "FontSize", 18);
text(ax1, 0.04, -0.048, '(a)', 'Units', 'normalized', 'FontSize', 17, 'FontWeight', 'bold');

% Plot observed and simulated values as lines
ax2 = axes('Position', [0.35, 0.15, 0.62, 0.75]); 
plot(uniqueMonths, monthlyQobs, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
hold on;
plot(uniqueMonths, monthlyQsim, '--r', 'LineWidth', 1.5);

% Add vertical reference lines
%xline(datetime(1995,1,1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
xline(datetime(1997,7,1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
%xline(datetime(2008,1,1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
xline(datetime(2010,7,1), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(ax2, 0.55, 0.7, 'Millennium Drought', 'Units', 'normalized', 'FontSize', 16, 'HorizontalAlignment', 'center');
%Add arrows pointing to the text annotations
ax2.FontSize = 14;
annotation('doublearrow', [0.471, 0.917], [0.65, 0.65], 'Color', 'k', 'LineWidth', 0.5);
text(ax2, -0.065, -0.15, '(b)', 'Units', 'normalized', 'FontSize', 17, 'FontWeight', 'bold');

% Customize plot appearance
xlim([datetime(1994,1,1), datetime(2011,12,31)]);
ax=gca;
set(ax, 'XTICK', datetime(1994:1:2011,1,1));
xtickformat(ax, 'yyyy');
xlabel('Year', 'FontSize', 17);
ylabel('Q (mm/month)', 'FontSize', 17);
legend("Observed streamflow", "Simulated streamflow by model Australia", 'Fontsize', 15, 'Location', 'northwest');
%grid on;
title('Monthly streamflow performance in catchment 406213', 'FontSize', 18, 'FontWeight', 'bold')
hold off;

if fileSave
    exportgraphics(gcf, [exp_path 'FigureS4.png']);
end
%% Figure S7
% Catchments to plot
catchments = {'406213', '406214', '407230'};
labels = {'(a)', '(b)', '(c)'};

% Make a wide white-background figure
figure;
set(gcf, 'Position', [100 200 1200 380]);

% Create a tiled layout with compact spacing
t = tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

for i = 1:length(catchments)
    ThisCatch = catchments{i};
    load(['Trade-offs/' ThisCatch, '_Tradeoffs.mat'])

    x = TradeoffsTable.equal_of3;
    y = TradeoffsTable.equal_Q;

    nexttile; hold on;
    scatter(x, y, [], TradeoffsTable.ed_equal, 'filled', 'HandleVisibility', 'off')
    clim([0, 1])  % consistent color scale for all panels

    xlim([0, 1])
    ylim([-1, 1])
    
    xlabel('Multi-annual storage dynamics (OF3)')
    if i == 1
        ylabel('Q performance (Mean of OF1 & OF2)')
    end
    text(0.04, -0.7, labels{i}, 'FontSize', 16, 'FontWeight', 'bold');

    grid on
    set(gca, 'FontSize', 13)
    title([ThisCatch ' equal weights'], 'FontSize', 17)
end

% Add a single colorbar spanning the whole layout
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'Distance to CTPP';
c.FontSize = 13;
cbarrow('up')

if fileSave
    exportgraphics(gcf, [exp_path 'FigureS7.png']);
end

%% Figure S8
% Model counts
with_count = 9;
without_count = 37;
start_angle = 90;
angle_with = 360 * with_count / (with_count + without_count);
angle_without = 360 - angle_with;

% Test data: [with_hist, with_syn, without_hist, without_syn]
test_data = [
    6, 2, 4, 1;
    6, 3, 1, 0;
    7, 4, 3, 1
];

% Colors
blue = [201 223 229]/255;
grey = [210 210 210]/255;

% Shifts
x_shift = 0.02;
y_shift = -0.02;

% Figure
f = figure;
f.Position = [600 300 1000 520];
annotation('line', [0.235 0.235] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
annotation('line', [0.475 0.475] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
annotation('line', [0.727 0.727] + x_shift, [0.2 0.9] + y_shift, ...
           'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);

annotation('textbox', [0.257, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '406213', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.505, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '406214', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.755, 0.62, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '407230', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');
annotation('textbox', [0.015, 0.22, 0.2, 0.05] + [x_shift y_shift 0 0], ...
    'String', '(HSC = hypothesized structual components)', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', ...
    'FontSize', 12);

%-----------------------No. models----------------------
annotation('textbox', [0, 0.68, 0.2, 0.05], ...
    'String', '9', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', [0.4, 0.6, 1], 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.05, 0.5, 0.2, 0.05], ...
    'String', '37', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', [130 130 130]/255, 'FontWeight','bold',...
    'FontSize', 15);

% 406213
annotation('textbox', [0.24, 0.85, 0.2, 0.05], ...
    'String', num2str(test_data(1,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.333, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(1,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.176, 0.27-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(1,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.23, 0.44-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(1,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

% 406214
annotation('textbox', [0.5, 0.858, 0.2, 0.05], ...
    'String', num2str(test_data(2,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.582, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(2,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.412, 0.385-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(2,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.48, 0.44-0.04, 0.2, 0.05], ...
    'String', num2str(test_data(2,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);
%407230
annotation('textbox', [0.727, 0.827, 0.2, 0.05], ...
    'String', num2str(test_data(3,1)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.833, 0.78, 0.2, 0.05], ...
    'String', num2str(test_data(3,2)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

annotation('textbox', [0.73, 0.395, 0.2, 0.05], ...
    'String', num2str(test_data(3,3)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontWeight','bold',...
    'FontSize', 15);
annotation('textbox', [0.665, 0.29, 0.2, 0.05], ...
    'String', num2str(test_data(3,4)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'Color', 'r', 'FontWeight','bold',...
    'FontSize', 15);

% -------- Combined pie subplot (left-middle) --------
combined_ax = axes('Position', [0.015 + x_shift, 0.4 + y_shift, 0.2, 0.45]);
hold(combined_ax, 'on'); axis(combined_ax, 'equal'); axis(combined_ax, 'off');

offset = 1.5; % for visual separation
fill_arc(0, 0, 1, start_angle + angle_with, start_angle + 360, grey, 'none'); % grey base
fill_arc(-0.07, 0.07, 1, start_angle - offset, start_angle + angle_with - offset, blue, 'none'); % blue offset

p1 = patch(NaN, NaN, [201 223 229]/255, 'EdgeColor', 'none'); % blue
p2 = patch(NaN, NaN, [210 210 210]/255, 'EdgeColor', 'none'); % grey

% Add legend
legend([p1 p2], {'include all HSC', 'not include all HSC'}, ...
    'Box', 'off', ...
    'Location', [0.055 + x_shift+0.01, 0.3 + y_shift, 0.1, 0.1], ...
    'FontSize', 13, 'FontWeight','bold');

% -------- Main 3-column layout --------
full_width = 0.2;
gap = 0.05;
x_start = 0.35 + x_shift;

for i = 1:3
    d = test_data(i, :);
    wh = d(1); ws = d(2); wth = d(3); wts = d(4);

    x_center = x_start + (i-1)*(full_width + gap);
    left_with = x_center - 0.2/4;
    left_without = x_center - 0.2/2;

    % --- Upper subplot (blue) ---
    ax1 = axes('Position', [left_with, 0.7 + y_shift, 0.2/2, 0.45/2]);
    hold(ax1, 'on'); axis(ax1, 'equal'); axis(ax1, 'off');
    fill_arc(0, 0, 1, start_angle, start_angle + angle_with, blue, 'none');
    fill_arc(0, 0, 1.05, start_angle, start_angle + angle_with * (wh / with_count), 'none', 'k');
    fill_arc(0, 0, 1.1, start_angle, start_angle + angle_with * (ws / with_count), 'none', 'r');

    % --- Lower subplot (grey) ---
    ax2 = axes('Position', [left_without, 0.19 + y_shift-0.04, 0.2, 0.45]);
    hold(ax2, 'on'); axis(ax2, 'equal'); axis(ax2, 'off');
    fill_arc(0, 0, 1, start_angle + angle_with, 360 + start_angle, grey, 'none');
    fill_arc(0, 0, 1.05, start_angle + angle_with, start_angle + angle_with + angle_without * (wth / without_count), 'none', 'k');
    if i == 3
        fill_arc(0, 0, 1.1, 1.3+start_angle + angle_with + angle_without * (wth / without_count), ...
            start_angle + angle_with + angle_without * (wth / without_count) + angle_without * (wts / without_count), 'none', 'r');
    else
        fill_arc(0, 0, 1.1, start_angle + angle_with, start_angle + angle_with + angle_without * (wts / without_count), 'none', 'r');
    end
end

% Create square patches with no fill and only colored edges
p1 = patch(NaN, NaN, 'w', 'EdgeColor', 'k', 'LineWidth', 1.2); % black edge only
p_spacer = patch(NaN, NaN, 'w', 'FaceColor', 'none', 'EdgeColor', 'none', 'Visible', 'off'); % invisible spacer
p2 = patch(NaN, NaN, 'w', 'EdgeColor', 'r', 'LineWidth', 1.2); % red edge only

% Add legend with two rows
legend([p1 p_spacer p2], ...
       {'Passed historic test', '', 'Passed synthetic test'}, ...
       'Box', 'off', ...
       'FontSize', 14, ...
       'FontWeight','bold',...
       'NumColumns', 3, ...
       'Position', [0.35 + x_shift, 0.12 + y_shift-0.02, 0.5, 0.03]);

if fileSave
    exportgraphics(gcf, [exp_path 'FigureS8.png']);
end


% ---------- Helper Function ----------
function fill_arc(xc, yc, r, theta1, theta2, facecolor, edgecolor)
    t = linspace(deg2rad(theta1), deg2rad(theta2), 100);
    x = [xc + r * cos(t), xc];
    y = [yc + r * sin(t), yc];
    if ischar(facecolor) && strcmp(facecolor, 'none')
        fill(x, y, 'w', 'FaceColor', 'none', 'EdgeColor', edgecolor, 'LineWidth', 1.2);
    else
        fill(x, y, facecolor, 'EdgeColor', edgecolor, 'LineWidth', 1.2);
    end
end

