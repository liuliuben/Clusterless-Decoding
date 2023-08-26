close all
%%
% Create a new figure and set its size
figure('Position', [100, 100, 1200, 600]);

% Plot rat's trajectory in black
plot(1:T, x, 'k-', 'LineWidth', 2);  % Increased line width

hold on;

% Assuming you have other plots that use colors
colors = parula(numCells);

% Improved axes appearance
axis tight;  % Tighten the axis
grid on;     % Display grid
set(gca, 'FontSize', 16, 'LineWidth', 1.5);  % Increase font size and axis line width
box on;  % Display a box around the plot

% Labels and title
set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
xlabel('Time (ms)', 'FontSize', 18,'FontWeight', 'bold');
ylabel('Rat Position', 'FontSize', 18,'FontWeight', 'bold');
title('Rat Position Over Time', 'FontSize', 20, 'FontWeight', 'bold');

% If you have other plots, you might want to introduce a legend
% legend({'Rat Trajectory', 'Other data'}, 'Location', 'best', 'FontSize', 14);

% Optional: Save the figure
% saveas(gcf, 'RatPositionOverTime.png');

%% jmi plotter

dm = .1; 
ms = 0:dm:15;
jmi = zeros(length(xs), length(ms));

for idx = 1:length(all_spike_times)
    jmi = jmi + (normpdf(xs', x(all_spike_times(idx)), sigma) * normpdf(ms, all_spike_marks(idx), sigm));
end

% --- First Plot ---
figure('Position', [100, 100, 1200, 600]);  % Define figure size
imagesc(xs, ms, jmi');
colormap('parula');  % Use jet colormap
colorbar;  % Add colorbar to indicate scale
set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
axis([-3 3 0 15]);
xlabel('Position', 'FontSize', 18, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
ylabel('Mark Value', 'FontSize', 18, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
title('Joint Mark Intensity', 'FontSize', 20, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here

% --- Second Plot ---

jointMarkIntensitiesPerCell = cell(1, numCells);

for i = 1:numCells
    cellSpikeTimes = spike_times_cell{i};
    cellSpikeMarks = spike_marks_cell{i};
    
    intensity = zeros(length(xs), length(ms));
    for idx = 1:length(cellSpikeTimes)
        intensity = intensity + (normpdf(xs', x(cellSpikeTimes(idx)), sigma) * normpdf(ms, cellSpikeMarks(idx), sigm));
    end
    jointMarkIntensitiesPerCell{i} = intensity / length(cellSpikeTimes);  % Normalize by the number of spikes
end

figure('Position', [100, 100, 1200, 600]);  % Define figure size
hold on;

% Plotting the first 10 cells with blue color
for i = 1:(numCells - 2)  % Assuming you have 12 cells and the last two are high amplitude
    intensity = jointMarkIntensitiesPerCell{i};
    placeField = mean(intensity, 2);
    plot(xs, placeField, 'b-', 'LineWidth', 1.5);  
end

% Plotting the last two cells with red color
for i = (numCells - 1):numCells  % For the last two cells
    intensity = jointMarkIntensitiesPerCell{i};
    placeField = mean(intensity, 2);
    plot(xs, placeField, 'r-', 'LineWidth', 2.5);
end

% Adjusting the plot
xlabel('Position (x)', 'FontSize', 18, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
ylabel('Intensity (Approx. Firing Rate)', 'FontSize', 18, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
title('Place Fields for Each Cell', 'FontSize', 20, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here
grid on;  % Add a grid
box on;  % Add a box around the plot

% Setting up the legend
hLow = plot(nan, nan, 'b-', 'LineWidth', 1.5);  % Dummy plot for legend
hHigh = plot(nan, nan, 'r-', 'LineWidth', 2.5);  % Dummy plot for legend
legend([hLow, hHigh], 'Low Amplitude Cells', 'High Amplitude Cells', 'FontSize', 14, 'Location', 'best');  % Adjusting the legend

hold off;

%% place fields example plot

dx = 0.01;  % Small step size for smooth curves
xs = -5:dx:5;  % Extended range for visualization

% Initialize a figure for place fields and set its size
figure('Position', [100, 100, 1200, 600]);
hold on;

% Plot properties
title('Place Fields for Place Cells', 'FontSize', 20,'FontWeight', 'bold');
xlabel('Position', 'FontSize', 18,'FontWeight', 'bold');
ylabel('Firing Rate', 'FontSize', 18,'FontWeight', 'bold');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');  % Add 'FontWeight', 'bold' here

% Improve axes and grid appearance
axis tight;  % Adjust the axes to fit the data
grid on;     % Display grid
set(gca, 'FontSize', 16, 'LineWidth', 1.5);  % Increase font size and axis line width
box on;  % Display a box around the plot
% ... [previous code remains unchanged]

% Data plotting
lowAmplitudePlaceFields = zeros(numLowAmplitudeCells, length(xs));

% Create empty arrays to hold plot handles
hLow = [];
hHigh = [];

for i = 1:numCells
    placeField = normpdf(xs, mu_x(i), sigma_x_values(i));
    
    % If you want to normalize
    % placeField = placeField / max(placeField);
    
    if i <= numLowAmplitudeCells
        lowAmplitudePlaceFields(i, :) = placeField;
        hTemp = plot(xs, placeField, 'b-', 'LineWidth', 1.5);  % Use a slightly bolder line for clarity
        hLow = [hLow; hTemp(1)];  % Store the first handle, if there's more than one
    else
        hTemp = plot(xs, placeField, 'r-', 'LineWidth', 2);
        hHigh = [hHigh; hTemp(1)];  % Store the first handle, if there's more than one
    end
end

% Calculate and plot the mean place field for low amplitude cells
meanLowAmplitudePlaceField = mean(lowAmplitudePlaceFields, 1);
hMean = plot(xs, meanLowAmplitudePlaceField, 'g-', 'LineWidth', 3);  % Even bolder for emphasis

% Add a legend using handles to ensure correct color association
legend([hLow(1), hHigh(1), hMean], 'Low Amplitude Cells', 'High Amplitude Cells', 'Mean Low Amplitude Place Field', 'FontSize', 14, 'Location', 'best');

hold off;

%%

% Set up a new figure with a larger default size
figure('Position', [100, 100, 1200, 800]);

% --- First Subplot (Rat Position Over Time with Spike Events) ---
subplot(2, 1, 1);
plot(1:T, x, 'k-', 'LineWidth', 2);  % Plot rat's trajectory in black with thicker line
hold on;

% Colors for different cells
colors = parula(numCells);

% Loop through each cell and plot its spikes
for cellIdx = 1:numCells
    cell_spike_times = spike_times_cell{cellIdx};
    
    scatter(cell_spike_times, x(cell_spike_times), 50, colors(cellIdx,:), 'filled', 'DisplayName', sprintf('Cell %d', cellIdx));
end

% Bold axes, labels, and title
xlabel('Time (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Rat Position', 'FontSize', 18, 'FontWeight', 'bold');
title('Rat Position Over Time with Spike Events', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');
legend('Location', 'northeastoutside');  % Place legend outside the plot to avoid overlap

% --- Second Subplot (Spike Amplitude Over Time) ---
subplot(2, 1, 2);
hold on;

for cellIdx = 1:numCells
    cell_spike_times = spike_times_cell{cellIdx};
    cell_spike_marks = spike_marks_cell{cellIdx};
    
    scatter(cell_spike_times, cell_spike_marks, 50, colors(cellIdx,:), 'filled', 'DisplayName', sprintf('Cell %d', cellIdx));
end

% Bold axes, labels, and title
xlabel('Spike Times (ms)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Spike Amplitude', 'FontSize', 18, 'FontWeight', 'bold');
title('Spike Amplitude Over Time', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');
legend('Location', 'northeastoutside');  % Place legend outside the plot to avoid overlap

% Tight layout to make use of all space

%% plot stats
% Set up figure parameters
figureSize = [100, 100, 1200, 600];  % Custom size for the figure
fontSize = 18;
lineWidth = 1.5;

% --- First Plot (Posterior Heatmap with Actual and Decoded Positions) ---
figure('Position', figureSize);
imagesc(1:T, xs, post);
set(gca, 'YDir', 'normal');  % Invert the y-axis to 'normal'
hold on;

% Plot Actual and Decoded Position
actualHandle = plot(1:T, x, 'w-', 'LineWidth', lineWidth);
decodedHandle = plot(1:T, decoded_x, 'k:', 'LineWidth', lineWidth);

xlabel('Time (ms)', 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel('Position', 'FontSize', fontSize, 'FontWeight', 'bold');
title('Posterior Heatmap with Actual and Decoded Positions', 'FontSize', fontSize + 2, 'FontWeight', 'bold');

% Manually create the legend
legend([actualHandle, decodedHandle], {'Actual Position', 'Decoded Position'}, 'Location', 'northeast');

set(gca, 'FontSize', fontSize, 'LineWidth', lineWidth, 'FontWeight', 'bold');
colorbar;
colormap('parula');

% --- Second Plot (Posterior Heatmap with Actual, Decoded Positions, and 95% HPD) ---
figure('Position', figureSize);
imagesc(1:T, xs, post);
set(gca, 'YDir', 'normal');  % Invert the y-axis to 'normal'
hold on;
plot(1:T, x, 'w-', 'LineWidth', lineWidth, 'DisplayName', 'Actual Position');
plot(1:T, decoded_x, 'k:', 'LineWidth', lineWidth, 'DisplayName', 'Decoded Position');
fill([1:T, fliplr(1:T)], [hpd_lower, fliplr(hpd_upper)], 'c', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', '95% HPD');

xlabel('Time (ms)', 'FontSize', fontSize, 'FontWeight', 'bold');
ylabel('Position', 'FontSize', fontSize, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', fontSize - 2);
title('Posterior Heatmap with Actual, Decoded Positions, and 95% HPD', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
set(gca, 'FontSize', fontSize, 'LineWidth', lineWidth, 'FontWeight', 'bold');
colorbar;
colormap('parula');
%%

% Set up figure parameters
figureSize = [100, 100, 800, 800];  % Custom size for the figure
fontSize = 18;
lineWidth = 1.5;
barWidth = 0.6;  % for a thicker bar
barColor = [0.7, 0.7, 0.7];  % Setting a gray color for the bars

% Initialize the figure with custom size
figure('Position', figureSize);

% --- Plotting RMSE ---
% subplot(1, 2, 1);

% Plot bar with error bar
% bar(1, RMSE_mean, 'BarWidth', barWidth, 'FaceColor', [0.7, 0.7, 0.7]);  % Setting a gray color for the bar
% hold on;
% errorbar(1, RMSE_mean, RMSE_SE, 'k.', 'LineWidth', lineWidth);
% 
% % Set labels, title, and adjust font sizes
% title('Root Mean Square Error', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
% ylabel('RMSE value', 'FontSize', fontSize, 'FontWeight', 'bold');
% xlabel('Trial', 'FontSize', fontSize, 'FontWeight', 'bold');
% ylim([0 1])
% set(gca, 'FontSize', fontSize, 'LineWidth', lineWidth, 'FontWeight', 'bold');

% --- Plotting mean HPDR width ---
subplot(1, 2, 1);

% Plot bar with error bar
bar(1, HPDR_width, 'BarWidth', barWidth, 'FaceColor', [0.7, 0.7, 0.7]);  % Setting a gray color for the bar
hold on;
errorbar(1, HPDR_width, width_SE, 'k.', 'LineWidth', lineWidth);

% Set labels, title, and adjust font sizes
title('Mean HPDR width', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
ylabel('x position', 'FontSize', fontSize, 'FontWeight', 'bold');
xlabel('Trial', 'FontSize', fontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', fontSize, 'LineWidth', lineWidth, 'FontWeight', 'bold');

% Plot bar with error bar
subplot(1, 2, 2);
bar(1, HPDR_mean, 'BarWidth', barWidth, 'FaceColor', barColor);
hold on;

errorbar(1, HPDR_mean, HPDR_SE, 'k.', 'LineWidth', lineWidth);
% Set labels, title, and adjust font sizes
title('HPDR across Trials', 'FontSize', fontSize + 2, 'FontWeight', 'bold');
ylabel('HPDR Value', 'FontSize', fontSize, 'FontWeight', 'bold');
xlabel('Metrics', 'FontSize', fontSize, 'FontWeight', 'bold');
set(gca, 'FontSize', fontSize, 'LineWidth', lineWidth, 'FontWeight', 'bold');
%%
dm = .1; 
ms = 0:dm:15;

% Assuming all_spike_times and all_spike_marks are cell arrays
num_cells = 12;

for c = 1:num_cells
    
    % Initialize a joint mark intensity for this cell
    jmi = zeros(length(xs), length(ms));
    
    % Extract the spike times and marks for this cell
    stc = spike_times_cell{c};
    stm = spike_marks_cell{c};
    
    for idx = 1:length(stc)
        jmi = jmi + (normpdf(xs', x(stc(idx)), sigma) * normpdf(ms, stm(idx), sigma));
    end

    % --- Plot for this cell ---
    figure('Position', [100, 100, 1200, 600]);  % Define figure size
    imagesc(xs, ms, jmi');
    colormap('parula');
    colorbar;
    set(gca, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1.5, 'FontWeight', 'bold');
    axis([-3 3 0 15]);
    xlabel('Position', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Mark Value', 'FontSize', 18, 'FontWeight', 'bold');
    title(['Joint Mark Intensity for Cell ' num2str(c)], 'FontSize', 20, 'FontWeight', 'bold');
    
    % Optionally add a pause or save each figure, if desired
    % pause(1); % Pauses for 1 second between plots
    % saveas(gcf, sprintf('Cell_%d.png', c)); % Saves each plot as a PNG
end

%%
figure;
hold on;
colors = lines(numCells);  % Different colors for different neurons

for cellIdx = 1:numCells
    plot(m_values, lambda_values(1000, :, cellIdx), 'Color', colors(cellIdx, :));  % Assuming time 500 as representative time point
end
xlabel('Mark Value');
ylabel('Firing Intensity');
title('Mark Field Representation');
