% -----------------------------------------------------
% Create safe snapshot matrix Y (handles NaT, NaN, Inf)
% Assumes Data exists
% -----------------------------------------------------

% Sort by time first
Data = sortrows(Data, 'datetime');

% Remove rows with invalid datetime or temperature
validRows = ~(isnat(Data.datetime) | isnan(Data.temperature) | isinf(Data.temperature));
Data = Data(validRows, :);

% Create 1-minute global time grid
% ---------------------------------------------------------
% Create a time grid of approximately 53,884 snapshots
% (≈ 53,884 minutes)
% Assumes Data.datetime exists and is sorted
% ---------------------------------------------------------

target_len = 500;    % desired number of snapshots

% Pick a valid random starting index
% ensure enough data is available after start_idx
start_idx = randi([1, height(Data) - target_len]);

% Define start and end times
t_start = Data.datetime(start_idx);
t_end   = t_start + minutes(target_len - 1);

% Create uniform 1-minute time grid
t_grid = (t_start : minutes(1) : t_end)';

m = length(t_grid);   % number of snapshots

disp(m);


% Initialize Y
Y = NaN(54, m);

for mote = 1:54

    % Extract data for this mote
    idx = (Data.moteid == mote);
    moteData = Data(idx, :);

    % Skip if no data for this mote
    if height(moteData) < 2
        warning('Mote %d has insufficient data. Filling with nearest values later.', mote);
        continue;
    end

    % Remove invalid timestamps or temperatures for this mote
    valid = ~(isnat(moteData.datetime) | isnan(moteData.temperature) | isinf(moteData.temperature));
    t_raw = moteData.datetime(valid);
    temp_raw = moteData.temperature(valid);

    % If still insufficient, skip
    if length(t_raw) < 2
        continue;
    end

    % Perform safe interpolation
    temp_interp = interp1(t_raw, temp_raw, t_grid, 'linear', 'extrap');

    % Store in snapshot matrix
    Y(mote, :) = temp_interp;
end

% Fill remaining NaNs using nearest values
Y = fillmissing(Y, 'nearest');

disp('Y created successfully:');
disp(size(Y));
clear idx m mot moteData t_end t_grid T_raw t_start temp_interp...
    temp_raw valid validRows t_raw mote

%%
sel = [Data.epoch,Data.moteid,Data.temperature];
moteId = Data.moteid;
sel2 = [];
for i = 1:length(sel)
    if sel(i,2)  < 55
        sel2 = sel(i,:);
    end
end
%%
sel1 = 
Range = round(linspace(1,65535,100));
a =[];
b = zeros(54,1);

a = find(sel(:,1) == Range(1));

