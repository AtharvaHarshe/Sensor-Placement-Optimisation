img = imread('lab.png');
figure; imshow(img); hold on;

% Load pixel coordinates (54×2)
% pixel_coords = [xpix, ypix];     % If already in workspace
% Otherwise: load pixel_coords.mat

% Sensor set to highlight
k = 21;
S = Pareto_Sets{k};     % e.g., [1 5 7 11 ...]

% Circle radius (adjust as needed)
radius = 22;
th = linspace(0, 2*pi, 200);

for i = 1:length(S)
    s = S(i);               % sensor index
    xc = pixel_coords(s,1) + radius*cos(th);
    yc = pixel_coords(s,2) + radius*sin(th);
    
    plot(xc, yc, 'r-', 'LineWidth', 3);
end

title(sprintf('Selected Sensors (k = %d) Highlighted', k), 'FontSize', 14);
hold off;
