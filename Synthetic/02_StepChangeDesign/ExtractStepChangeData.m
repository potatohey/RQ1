load("HistoricData.mat");

Year = 2010;
Catch = 'x407230';
wetYears = 10; % how many wet years per cycle
dryYears = 10; % how many dry years per cycle
cycle = 2; % 2*(10+10) = 40 years in total
wet2dry = [1 1 1 0.5]; % shift 50% down

P_slice = HistoricData.precip_daily(HistoricData.precip_daily.year == Year, {'year', 'month','day', Catch});
PET_slice = HistoricData.PET_daily(HistoricData.PET_daily.year == Year, {'year', 'month','day', Catch});

P_wet = repmat(P_slice, wetYears, 1);
P_dry = repmat(P_slice.*wet2dry, dryYears, 1);

P = [P_wet; P_dry];
PET = repmat(PET_slice, wetYears+dryYears, 1);

P = repmat(P, cycle, 1);
PET = repmat(PET, cycle, 1);

P.year = ones(size(P, 1), 1);

for i = 2:size(P, 1)
    if P{i, 'month'} == 1 && P{i, 'day'} == 1
        P{i, 'year'} = P{i-1, 'year'} + 1;
    else
        P{i, 'year'} = P{i-1, 'year'};
    end
end

PET.year = P.year;

StepChangeData = struct("P", P, "PET", PET);
StepChangeData.P.Properties.VariableNames(1) = "Year";
StepChangeData.P.Properties.VariableNames(2) = "Month";
StepChangeData.P.Properties.VariableNames(3) = "Day";
StepChangeData.P.Properties.VariableNames(4) = "P";
StepChangeData.PET.Properties.VariableNames(1) = "Year";
StepChangeData.PET.Properties.VariableNames(2) = "Month";
StepChangeData.PET.Properties.VariableNames(3) = "Day";
StepChangeData.PET.Properties.VariableNames(4) = "PET";


