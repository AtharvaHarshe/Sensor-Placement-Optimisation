sel = [Data.epoch , Data.moteid , Data.temperature];

% Find all unique sensors (should be 54)
sensorList = unique(sel(:,2));

% Create a cell array to store clusters
clusters = {};

for i = 1:length(sensorList)
    sensorID = sensorList(i);
  
    % Extract rows belonging to this sensor
    rows = sel(sel(:,2) == sensorID, :);
    
    % Store epoch + temp only
    clusters{i,1} = rows(:, [1 3]);   % [epoch, temperature]
  
end