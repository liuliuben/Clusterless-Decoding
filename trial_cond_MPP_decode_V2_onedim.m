clear
close all
%% cmon 
numTrials = 10;

% Define joint mark intensity function
lambda_t_m = @(x_pos, m, mu_x_cell, sigma_x_cell, mu_m_cell, sigma_m_cell, inv_matrix, alphaValue_cell) ...
    alphaValue_cell * exp(-0.5 * ([x_pos - mu_x_cell; m - mu_m_cell]' ...
    * inv_matrix ...
    * [x_pos - mu_x_cell; m - mu_m_cell]));

% Preallocate arrays for results
RMSE_values_condition1 = zeros(1, numTrials);
Width_condition1 = zeros(1, numTrials);
HPDR_values_condition1 = zeros(1, numTrials);
RMSE_values_condition2 = zeros(1, numTrials);
Width_condition2 = zeros(1,numTrials);
HPDR_values_condition2 = zeros(1,numTrials);
RMSE_values_condition3 = zeros(1,numTrials);
Width_condition3 = zeros(1,numTrials);
HPDR_values_condition3 = zeros(1,numTrials);

for trial = 1:numTrials
    
    stored_high_amplitude_spike_times = [];
    stored_high_amplitude_spike_marks = [];
    
    % Simulate rat's movement
    T = 1000; % Total time
    alpha = 0.98; % Rate of decay or persistence of the rat's motion in the previous direction
    sigma = 0.15; % Standard deviation of the noise in the rat's motion
    x = zeros(1, T); % position vector
   
    % Simulate rat movement across a linear track using a 1D random walk model
    for t = 2:T
        x(t) = alpha * x(t-1) + sigma * randn();  
    end

    % Simulate spiking activity of neurons
    % Define properties for neurons with higher amplitude marks
    mu_x_high = [-1.5, 1.5];
    mu_m_high = [11, 13];
    
    for cond = 1:3
        
        tic
        lowcells = [10,0,0];
        numLowAmplitudeCells = lowcells(cond);
        numHighAmplitudeCells = 2;
        numCells = numLowAmplitudeCells + numHighAmplitudeCells;
        sigma_x_high = sqrt(0.1) * ones(1, numHighAmplitudeCells);
        sigk = sigma_x_high(1); %this is just for the clustered plotting
        sigma_m_values = 2 * ones(1, numCells);  
        alphaValues = 100 * ones(1, numCells);  
        
        if cond == 1 % low and high amplitude
            pos_range = [-3, 3]; % Randomly generate the centers (mu_x_low) for the place fields so they span the entire position grid
            mu_spacing = (pos_range(2) - pos_range(1)) / numLowAmplitudeCells;
            mu_x_low = linspace(pos_range(1) + mu_spacing/2, pos_range(2) - mu_spacing/2, numLowAmplitudeCells); 
            jitter_range = 0.2;  
            jitter_values = (rand(1, numLowAmplitudeCells) - 0.5) * 2 * jitter_range; % This will give values between -jitter_range to +jitter_range
            mu_x_low = mu_x_low + jitter_values;
            % Define the range of standard deviations
            sigma_range = [0.5, 1.5];  % Adjust these values for desired width range
            % Generate random standard deviations for each cell, such that they are wider
            sigma_x_low = sigma_range(1) + (sigma_range(2) - sigma_range(1)) * rand(1, numLowAmplitudeCells); 
            mu_m_low = 2 + rand(1, numLowAmplitudeCells);  % Random marks between 2 and 3
            
            % Combine properties
            mu_x = [mu_x_low, mu_x_high];
            mu_m = [mu_m_low, mu_m_high];
            sigma_x_values = [sigma_x_low, sigma_x_high];
            sigm = sigma_m_values(1);
                                
            % Calculate joint intensity function values 
            % The lambda_values array stores the joint intensity of spiking for every combination of time, mark value, and neuron.
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

            % Simulate spiking activity for each neuron
            % For each neuron and each time point, a random number is compared to the total intensity 
            % (summed across all possible mark values). If the random number is smaller, a spike occurs. 
            % The mark value of the spike is then determined probabilistically based on the joint intensity values.
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
            fprintf('Time taken for low amp initiation: %.4f seconds\n', toc);

            % Generate a global list of spike times and amplitudes across all cells.
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

            % Decode with JMI
            dx = 1e-2; xs = -4:dx:4;
            OSM = normpdf(xs'*ones(size(xs)),ones(size(xs'))*alpha*xs,sigma);
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
            
             % Storing high amplitude spike times and marks for later conditions
            stored_high_amplitude_spike_times = [spike_times_cell{end-1}, spike_times_cell{end}]; % Last two cells are high amplitude cells
            stored_high_amplitude_spike_marks = [spike_marks_cell{end-1}, spike_marks_cell{end}]; 

            [~, maxIdx] = max(post, [], 1);
            decoded_x = xs(maxIdx);
            error = x - decoded_x;
            RMSE = sqrt(mean(error.^2))

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
                                          
        elseif cond == 2 %use only high amplitude
            
            all_spike_times = stored_high_amplitude_spike_times;
            all_spike_marks = stored_high_amplitude_spike_marks;
            % Sort the spikes in time order.
            [all_spike_times, sortIdx] = sort(all_spike_times);
            all_spike_marks = all_spike_marks(sortIdx);

            % Decode with JMI
            dx = 1e-2; xs = -4:dx:4;
            OSM = normpdf(xs'*ones(size(xs)),ones(size(xs'))*alpha*xs,sigma);
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

            [~, maxIdx] = max(post, [], 1);
            decoded_x = xs(maxIdx);
            error = x - decoded_x;
            RMSE = sqrt(mean(error.^2))

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
               
        elseif cond == 3
            
            all_spike_times = stored_high_amplitude_spike_times;
            all_spike_marks = stored_high_amplitude_spike_marks;
            % Sort the spikes in time order.
            [all_spike_times, sortIdx] = sort(all_spike_times);
            all_spike_marks = all_spike_marks(sortIdx);
            
            id = kmeans(all_spike_marks',2); 
            id1 = (id==1)'; 
            id2 = (id==2)';
            spike_positions_id1 = x(all_spike_times(id1));
            spike_positions_id2 = x(all_spike_times(id2));
            
            % Calculate and plot lambda values for each cluster
            lambda1_high = sum(normpdf(xs'*ones(1, numel(spike_positions_id1)), ones(size(xs'))*spike_positions_id1, sigk), 2) ./ sum(normpdf(xs'*ones(size(x)), ones(size(xs'))*x, sigk), 2);
            lambda2_high = sum(normpdf(xs'*ones(1, numel(spike_positions_id2)), ones(size(xs'))*spike_positions_id2, sigk), 2) ./ sum(normpdf(xs'*ones(size(x)), ones(size(xs'))*x, sigk), 2);
            postc = normpdf(xs'*ones(1,T),0,1); %initialize probabilities           
            OSMc = normpdf(xs'*ones(size(xs)),ones(size(xs'))*alpha*xs,sigma);
            
            for t = 2:T
                os = OSMc * postc(:, t-1); 
                os = os / sum(os) / dx;

                % Check for spikes at time t for both clusters
                spike_at_t_id1 = any(all_spike_times(id1) == t);
                spike_at_t_id2 = any(all_spike_times(id2) == t);

                l = lambda1_high .^ spike_at_t_id1 .* exp(-lambda1_high) .* lambda2_high .^ spike_at_t_id2 .* exp(-lambda2_high);

                postc(:, t) = os .* l; 
                postc(:, t) = postc(:, t) / sum(postc(:, t)) / dx;
            end
            fprintf('Time taken for decoding: %.4f seconds\n', toc);

            [~, maxIdxc] = max(postc, [], 1);
            decoded_position = xs(maxIdxc);
            errorc = x - decoded_position;
            RMSEc = sqrt(mean(errorc.^2))
            
            % Initialize variables to store HPD intervals
            hpd_lowerc = zeros(1, T);
            hpd_upperc = zeros(1, T);
            coveragec = zeros(1, T);

            for t = 1:T
                % Sort the posterior for time step t and get indices
                [sorted_postc, sorted_indicesc] = sort(postc(:,t), 'descend');
                % Find the shortest interval with cumsum >= 0.95
                cum_sorted_postc = cumsum(sorted_postc) * dx;
                hpdc = find(cum_sorted_postc >= 0.95, 1);
                % Convert the HPD from indices to positions
                hpd_indicesc = sorted_indicesc(1:hpdc);
                hpd_lowerc(t) = min(xs(hpd_indicesc));
                hpd_upperc(t) = max(xs(hpd_indicesc));
                % Determine if true position is within the HPD
                coveragec(t) = x(t) >= hpd_lowerc(t) && x(t) <= hpd_upperc(t);
            end
            % Calculate the width of the HPD region for each time step
            hpd_widthsc = hpd_upperc - hpd_lowerc;
            % Calculate the average width
            average_hpd_widthc = mean(hpd_widthsc);
            % Calculate the percentage of time the true position is covered by the HPD
            percentage_coveragec = mean(coveragec) * 100;             
        end

        if cond == 1
            % Store results for condition 1
            RMSE_values_condition1(trial) = RMSE;
            Width_condition1(trial) = average_hpd_width;
            HPDR_values_condition1(trial) = percentage_coverage;
        elseif cond == 2
            % Store results for condition 2
            RMSE_values_condition2(trial) = RMSE;
            Width_condition2(trial) = average_hpd_width;
            HPDR_values_condition2(trial) = percentage_coverage;
        elseif cond == 3
            % Store results for condition 3
            RMSE_values_condition3(trial) = RMSEc;
            Width_condition3(trial) = average_hpd_widthc;
            HPDR_values_condition3(trial) = percentage_coveragec;
        end   
    end
end

% Compute mean and SE
% Condition 1
RMSE_mean_1 = mean(RMSE_values_condition1);
RMSE_SE_1 = std(RMSE_values_condition1) / sqrt(numTrials);
width_mean_1 = mean(Width_condition1);
width_SE_1 = std(Width_condition1) / sqrt(numTrials);
HPDR_mean_1 = mean(HPDR_values_condition1);
HPDR_SE_1 = std(HPDR_values_condition1) / sqrt(numTrials);

% Condition 2
RMSE_mean_2 = mean(RMSE_values_condition2);
RMSE_SE_2 = std(RMSE_values_condition2) / sqrt(numTrials);
width_mean_2 = mean(Width_condition2);
width_SE_2 = std(Width_condition2) / sqrt(numTrials);
HPDR_mean_2 = mean(HPDR_values_condition2);
HPDR_SE_2 = std(HPDR_values_condition2) / sqrt(numTrials);

% Condition 2
RMSE_mean_3 = mean(RMSE_values_condition3);
RMSE_SE_3 = std(RMSE_values_condition3) / sqrt(numTrials);
width_mean_3 = mean(Width_condition3);
width_SE_3 = std(Width_condition3) / sqrt(numTrials);
HPDR_mean_3 = mean(HPDR_values_condition3);
HPDR_SE_3 = std(HPDR_values_condition3) / sqrt(numTrials);

figure;

% RMSE
subplot(3, 1, 1);
bar([1 2 3], [RMSE_mean_1 RMSE_mean_2 RMSE_mean_3]);
hold on;
errorbar([1 2 3], [RMSE_mean_1 RMSE_mean_2 RMSE_mean_3], [RMSE_SE_1 RMSE_SE_2 RMSE_SE_3], 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Condition 1', 'Condition 2', 'Condition 3'});
title('RMSE Comparison');
ylabel('RMSE Value');

% Width
subplot(3, 1, 2);
bar([1 2 3], [width_mean_1 width_mean_2 width_mean_3]);
hold on;
errorbar([1 2 3], [width_mean_1 width_mean_2 width_mean_3], [width_SE_1 width_SE_2 width_SE_3], 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Condition 1', 'Condition 2', 'Condition 3'});
title('RMSE Comparison');
ylabel('RMSE Value');

% HPDR
subplot(3, 1, 3);
bar([1 2 3], [HPDR_mean_1 HPDR_mean_2 HPDR_mean_3]);
hold on;
errorbar([1 2 3], [HPDR_mean_1 HPDR_mean_2 HPDR_mean_3], [HPDR_SE_1 HPDR_SE_2 HPDR_SE_3], 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Condition 1', 'Condition 2', 'Condition 3'});
title('HPDR Comparison');
ylabel('HPDR Value');
