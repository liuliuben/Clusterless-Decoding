%%
%across 100 trials average high amp = 60 and low amp = 1800;
numTrials = 100;
store1 = [];
store2 = [];
% Define joint mark intensity function
lambda_t_m = @(x_pos, m, mu_x_cell, sigma_x_cell, mu_m_cell, sigma_m_cell, inv_matrix, alphaValue_cell) ...
    alphaValue_cell * exp(-0.5 * ([x_pos - mu_x_cell; m - mu_m_cell]' ...
    * inv_matrix ...
    * [x_pos - mu_x_cell; m - mu_m_cell]));

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
    mu_x_high = [-1.5, 1.5];
    mu_m_high = [11, 13];
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
    
    num_spikes_lowamp_neurons = 0;
    for neuron_idx = 1:10
        num_spikes_lowamp_neurons = num_spikes_lowamp_neurons + length(spike_times_cell{neuron_idx});
    end

    num_spikes_highamp_neurons = length(spike_times_cell{end}) + length(spike_times_cell{end-1});
    store1(trial) = num_spikes_highamp_neurons;
    store2(trial) = num_spikes_lowamp_neurons;
end 
mean(store1)
mean(store2)

% Calculate the standard deviation
std_high_amp = std(store1);
std_low_amp = std(store2);

% Calculate the SEM
sem_high_amp = std_high_amp / sqrt(numTrials);
sem_low_amp = std_low_amp / sqrt(numTrials);

mean_high_amp = mean(store1);
mean_low_amp = mean(store2);

figure;
bar([1,2], [mean_high_amp, mean_low_amp]);
hold on;
errorbar([1,2], [mean_high_amp, mean_low_amp], [sem_high_amp, sem_low_amp], '.');
xticks([1,2]);
xticklabels({'High Amp', 'Low Amp'});
ylabel('Spike Count');
title('Average Spike Counts with Error Bars');
