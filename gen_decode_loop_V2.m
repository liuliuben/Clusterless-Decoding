clear;
% close all;

% Define joint mark intensity function
    lambda_t_m = @(x_pos, m, mu_x_cell, sigma_x_cell, mu_m_cell, sigma_m_cell, inv_matrix, alphaValue_cell) ...
        alphaValue_cell * exp(-0.5 * ([x_pos - mu_x_cell; m - mu_m_cell]' ...
        * inv_matrix ...
        * [x_pos - mu_x_cell; m - mu_m_cell]));

numTrials = 1;
RMSE_values = [];
HPDR_values = []; 

for trial = 1:numTrials
    
    alpha = 0.98; % Rate of decay or persistence of the rat's motion in the previous direction
    sigma = 0.15; % Standard deviation of the noise in the rat's motion
    T = 1000; % Total time
    x = zeros(1, T); % position vector

    % Simulate rat movement across a linear track using a 1D random walk model
    for t = 2:T
        x(t) = alpha * x(t-1) + sigma * randn();  
    end

    tic;
    % Define number of cells
    numLowAmplitudeCells = 10;
    numHighAmplitudeCells = 2;
    numCells = numLowAmplitudeCells + numHighAmplitudeCells;

    % Define properties for neurons with higher amplitude marks
    mu_x_high = [-1.2, 1.2];
    mu_m_high = [11, 13]; %close to 11 and 13 - 10/12
    sigma_x_high = sqrt(0.1) * ones(1, numHighAmplitudeCells);
    sigk = sigma_x_high(1); %this is just for the clustered plotting
    % Randomly generate properties for low amplitude cells
    % Randomly generate the centers (mu_x_low) for the place fields so they span the entire position grid
    pos_range = [-3, 3]; % Assuming the position grid ranges from -3 to 3
    mu_spacing = (pos_range(2) - pos_range(1)) / numLowAmplitudeCells;
    mu_x_low = linspace(pos_range(1) + mu_spacing/2, pos_range(2) - mu_spacing/2, numLowAmplitudeCells); 
    jitter_range = 0.2;  
    jitter_values = (rand(1, numLowAmplitudeCells) - 0.5) * 2 * jitter_range; % This will give values between -jitter_range to +jitter_range
    mu_x_low = mu_x_low + jitter_values;
    % Define the range of standard deviations
    sigma_range = [0.5, 1.5];  % Adjust these values for desired width range
    % Generate random standard deviations for each cell, such that they are wider
    sigma_x_low = sigma_range(1) + (sigma_range(2) - sigma_range(1)) * rand(1, numLowAmplitudeCells); 

%     mu_x_low = rand(1, numLowAmplitudeCells) * 4 - 2  % Random place fields between -2 and 2
    mu_m_low = 2 + rand(1, numLowAmplitudeCells);  % Random marks between 2 and 3
%     sigma_x_low = 0.5 + rand(1, numLowAmplitudeCells) * 1.5;  % Random width of each place field between 0.5 and 2
    
    % Combine properties
    mu_x = [mu_x_low, mu_x_high];
    mu_m = [mu_m_low, mu_m_high];
    sigma_x_values = [sigma_x_low, sigma_x_high];
    alphaValues = 100 * ones(1, numCells);  
    sigma_m_values = 2 * ones(1, numCells);  
    sigm = sigma_m_values(1);

    fprintf('Time taken for rat movement simulation: %.4f seconds\n', toc);
    
    tic;
    
    % Calculate joint intensity function values 
    % The lambda_values array stores the joint intensity of spiking for every combination of time, 
    % mark value, and neuron.
    markinc = 0.1;
    m_values = 0.1:markinc:15;

    lambda_values = zeros(T, length(m_values), numCells);
    for cellIdx = 1:numCells
        inv_matrix = inv([sigma_x_values(cellIdx)^2, 0; 0, sigma_m_values(cellIdx)^2]);
        for t = 1:T
            x_pos = x(t);  % Retrieve rat's position at time t
            for m_idx = 1:length(m_values)
                lambda_values(t, m_idx, cellIdx) = lambda_t_m(x_pos, m_values(m_idx), mu_x(cellIdx), sigma_x_values(cellIdx), mu_m(cellIdx), sigma_m_values(cellIdx), inv_matrix, alphaValues(cellIdx));
            end
        end
    end
    fprintf('Time taken for lambda calculations: %.4f seconds\n', toc);

    % Simulate spiking activity for each neuron
    % For each neuron and each time point, a random number is compared to the total intensity 
    % (summed across all possible mark values). If the random number is smaller, a spike occurs. 
    % The mark value of the spike is then determined probabilistically based on the joint intensity values.
    tic;
    
    spike_times_cell = cell(numCells, 1);
    spike_marks_cell = cell(numCells, 1);
    for cellIdx = 1:numCells
        spike_times = [];
        spike_marks = [];
        for t = 1:T
            rate = sum(lambda_values(t, :, cellIdx));
            if rand() < rate*0.0001
                spike_times = [spike_times, t];
                [~, idx] = max(rand(1, length(m_values)) < cumsum(lambda_values(t, :, cellIdx) / rate));
                spike_marks = [spike_marks, m_values(idx)];
            end
        end
        spike_times_cell{cellIdx} = spike_times;
        spike_marks_cell{cellIdx} = spike_marks;
    end

    if length(spike_times_cell) < 2
        error('The cell array contains less than two neurons.');
    end

    num_spikes_highamp_neurons = length(spike_times_cell{end}) + length(spike_times_cell{end-1})

    fprintf('Time taken for spike simulation: %.4f seconds\n', toc);
    
    %% Generate a global list of spike times and amplitudes across all cells.
    tic;
    all_spike_times = [];
    all_spike_marks = [];
    for cellIdx = 1:numCells
        all_spike_times = [all_spike_times, spike_times_cell{cellIdx}];
        all_spike_marks = [all_spike_marks, spike_marks_cell{cellIdx}];
    end
    % Sort the spikes in time order.
    [all_spike_times, sortIdx] = sort(all_spike_times);
    all_spike_marks = all_spike_marks(sortIdx);

    %% Decode with JMI
    dx = 1e-2; xs = -4:dx:4;
    post = normpdf(xs'*ones(1,T),0,1); 
    OSM = normpdf(xs'*ones(size(xs)),ones(size(xs'))*alpha*xs,sigma);

    %% Decode with JMI
    post = normpdf(xs'*ones(1,T),0,1); %initialize probabilities 
    % Compute the ground intensity
    groundIntensity = sum(normpdf(xs'*ones(size(all_spike_times)), ones(size(xs'))*x(all_spike_times), sigma),2)./sum(normpdf(xs'*ones(size(x)), ones(size(xs'))*x, sigma),2);

    for t=2:T 
        os = OSM*post(:,t-1); 
        os = os/sum(os)/dx; % step1 evolution matrix A
        l = exp(-groundIntensity); % step2

        % Find the spikes for the current time
        spikeIndicesForTimeT = find(all_spike_times == t);

        if ~isempty(spikeIndicesForTimeT)
            productOfSpikes = ones(length(xs), 1);  % Initializing to ones to ensure proper pointwise multiplication
            for spikeIdx = spikeIndicesForTimeT
                likelihoodPerSpike = (normpdf(xs'*ones(size(all_spike_times)), ones(size(xs'))*x(all_spike_times),sigma) * normpdf(all_spike_marks, all_spike_marks(spikeIdx), sigm)'); % multiply for the second dimension of marks * normpdf(all_spike_marks2, all_spike_marks(spikeIdx), sigm)'
                likelihoodPerSpike = likelihoodPerSpike ./ sum(normpdf(xs'*ones(size(x)),ones(size(xs'))*x,sigma),2);
                productOfSpikes = productOfSpikes .* likelihoodPerSpike;
            end
            l = l .* productOfSpikes;
        end
        post(:,t) = os.*l; post(:,t) = post(:,t)/sum(post(:,t))/dx;
    end
    fprintf('Time taken for decoding: %.4f seconds\n', toc);
    
    tic
    [~, maxIdx] = max(post, [], 1);
    decoded_x = xs(maxIdx);
    error = x - decoded_x;
    RMSE = sqrt(mean(error.^2));

    % Initialize variables to store HPD intervals
    hpd_lower = zeros(1, T);
    hpd_upper = zeros(1, T);
    coverage = zeros(1, T);

    for t = 1:T
        % Sort the posterior for time step t and get indices
        [sorted_post, sorted_indices] = sort(post(:,t), 'descend');

        % Find the shortest interval with cumsum >= 0.95
        cum_sorted_post = cumsum(sorted_post) * dx;
        hpd = find(cum_sorted_post >= 0.95, 1);

        % Convert the HPD from indices to positions
        hpd_indices = sorted_indices(1:hpd);
        hpd_lower(t) = min(xs(hpd_indices));
        hpd_upper(t) = max(xs(hpd_indices));

        % Determine if true position is within the HPD
        coverage(t) = x(t) >= hpd_lower(t) && x(t) <= hpd_upper(t);
    end
    % Calculate the width of the HPD region for each time step
    hpd_widths = hpd_upper - hpd_lower;

    % Calculate the average width
    average_hpd_width = mean(hpd_widths);
    % Calculate the percentage of time the true position is covered by the HPD
    percentage_coverage = mean(coverage) * 100;
 
    fprintf('Time taken for stats: %.4f seconds\n', toc);
    
    fprintf('Root Mean Square Error (RMSE): %.4f\n', RMSE);
    fprintf('hpdr average width: %.4f\n', average_hpd_width);
    fprintf('hpdr mean coverage: %.4f\n', percentage_coverage);

    RMSE_values(trial) = RMSE;
    HPDR_values(trial) = percentage_coverage;
    
    % ... [Your code up to this point remains unchanged]

% Compute the mean trajectory for low amplitude neurons
mean_low_amplitude_trajectory = mean(cell2mat(spike_times_cell(1:numLowAmplitudeCells)'),2);  % Assuming your data is in columns

% Compute the second derivative
second_derivative = diff(mean_low_amplitude_trajectory, 2);

% Find 100ms windows for highest and lowest second derivative
[~, idx_max] = max(second_derivative);
[~, idx_min] = min(second_derivative);
window_size = 100;  % 100ms window

% Make sure the indices are valid for the window size
if idx_max - window_size/2 < 1
    idx_max = window_size/2 + 1;
end
if idx_max + window_size/2 > T
    idx_max = T - window_size/2;
end
if idx_min - window_size/2 < 1
    idx_min = window_size/2 + 1;
end
if idx_min + window_size/2 > T
    idx_min = T - window_size/2;
end

% Extract the indices for the two windows
window_max = (idx_max - window_size/2):(idx_max + window_size/2);
window_min = (idx_min - window_size/2):(idx_min + window_size/2);

% Compute RMSE for each window
error_max = x(window_max) - decoded_x(window_max);
RMSE_max = sqrt(mean(error_max.^2));
error_min = x(window_min) - decoded_x(window_min);
RMSE_min = sqrt(mean(error_min.^2));

fprintf('Root Mean Square Error (RMSE) for high second derivative window: %.4f\n', RMSE_max);
fprintf('Root Mean Square Error (RMSE) for low second derivative window: %.4f\n', RMSE_min);

end

% Compute mean and SE
RMSE_mean = mean(RMSE_values);
RMSE_SE = std(RMSE_values) / sqrt(numTrials);
HPDR_width = average_hpd_width;
width_SE = std(HPDR_width)/ sqrt(numTrials);
HPDR_mean = mean(HPDR_values);
HPDR_SE = std(HPDR_values) / sqrt(numTrials);
fprintf('rmse trial width: %.4f\n', RMSE_mean);
fprintf('hpdr trial coverage: %.4f\n', HPDR_mean);

% Plotting
figure;
% Plotting RMSE
subplot(3, 1, 1);
bar(1, RMSE_mean);
ylim([0 1])
hold on;
errorbar(1, RMSE_mean, RMSE_SE, 'k.', 'LineWidth', 1.5);
title('RMSE across trials');
ylabel('RMSE value');

subplot(3, 1, 2);
bar(1, HPDR_width);
hold on;
errorbar(1, HPDR_width, width_SE, 'k.', 'LineWidth', 1.5);
title('mean HPDR width');
ylabel('width (x position)');

% Plotting HPDR
subplot(3, 1, 3);
bar(1, HPDR_mean);
hold on;
errorbar(1, HPDR_mean, HPDR_SE, 'k.', 'LineWidth', 1.5);
title('HPDR across trials');
ylabel('HPDR value');
xlabel('Metrics');