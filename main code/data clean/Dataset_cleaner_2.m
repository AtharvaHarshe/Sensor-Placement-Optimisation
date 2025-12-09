%clc; clear all; close all;

clusters2 = clusters;

for i = 1:length(clusters2)

    temp = clusters2{i};     % Nx2 matrix [epoch, temp]
    T = temp(:,2);           % temperature column
    
    % Logical filter for GOOD data
    good = (T >= 0) & (T <= 50) & (T ~= 122.153);
    
    % Apply filter → remove bad rows
    clusters2{i} = temp(good, :);
   
end
avg = [];
for i = 1:length(clusters2)

    temp = clusters2{i};     % Nx2 matrix [epoch, temp]
    avg(i,1) = mean(temp(:,2));           % temperature column
        
   
end
%%
% Number of sensors
numSensors = length(clusters2);

% Preallocate final result
snapshot = cell(numSensors, 1);

% The 100 uniformly spaced target epochs
Range = round(linspace(1, 2^16, 100));

for s = 1:numSensors
    
    temp = clusters2{s};        % Nx2 [epoch , temperature]
    avgValue = avg(s);    % average temperature for this sensor
    
    snap = zeros(100, 2);       % will store 100 rows [epoch , temp]
    
    for i = 1:100
        target = Range(i);
        
        % ---- 1. Check exact epoch match ----
        idx = find(temp(:,1) == target);
        
        % ---- 2. If not found, check ±1, ±2 ----
        if isempty(idx)
            offsets = [1 2 -1 -2 3 4 -3 -4 ];
            for k = 1:length(offsets)
                idx = find(temp(:,1) == (target + offsets(k)));
                if ~isempty(idx)
                    break;
                end
            end
        end
        
        % ---- 3. If still not found, assign average ----
        if isempty(idx)
            snap(i,:) = [target , avgValue];
        else
            snap(i,:) = temp(idx(1), :);   % use first found match
        end
    end
    
    snapshot{s} = snap;   % store the 100 readings for this sensor
end
snapshot = snapshot(1:54);
numSensors = length(snapshot);   % should be 54

TrainSet = zeros(numSensors, 70);
TestSet  = zeros(numSensors, 30);

for s = 1:numSensors
    
    tempData = snapshot{s}(:,2);   % extract temperature column (100x1)
    
    % First 70 → training
    TrainSet(s, :) = tempData(1:70);
    
    % Last 30 → testing
    TestSet(s, :) = tempData(71:100);
end


