function [stats, fh] = ephysio_focality_inertia(data, Y_vars,Y_labels, nperms)


if ischar(Y_vars); Y_vars = {Y_vars}; end
if ischar(Y_labels); Y_labels = {Y_labels}; end


% Prepare for plotting
fh = figure("Position",[100 100 numel(Y_vars)*300 300]);

x_mni = data.x_flip;
y_mni = data.y_flip;
z_mni = data.z_flip;

stats = table();

% Loop through each variable in Y_vars
for idx = 1:length(Y_vars)
    power_variable = Y_vars{idx};  % Current variable
    pow_label = Y_labels{idx};

    % Initialize the figure for each variable (create subplot for each)
    
    ns = numSubplots(numel(Y_vars));
    subplot(ns(1),ns(2), idx);


    % Dynamic weights assignment for current power variable
    weights = data.(power_variable);  % Dynamic column reference

    % Calculate weighted average of DBS coordinates
    x_avg = sum(x_mni .* weights) / sum(weights);
    y_avg = sum(y_mni .* weights) / sum(weights);
    z_avg = sum(z_mni .* weights) / sum(weights);
    weighted_mean = [x_avg, y_avg, z_avg];

    [max_val, max_idx] = max(weights);

    fprintf("Variable: %s\n", power_variable);
    fprintf("  Weighted Mean: [%.4f, %.4f, %.4f]\n", weighted_mean);
    fprintf("  Max Value:     %.4f (row %d)\n", max_val, max_idx);
    fprintf("  Coordinate at Max: [%.4f, %.4f, %.4f]\n", x_mni(max_idx),y_mni(max_idx),z_mni(max_idx));

    % Compute Euclidean distances to the weighted average
    euclidean_distances = sqrt((x_mni - x_avg).^2 + (y_mni- y_avg).^2 + (z_mni - z_avg).^2);

    % Calculate weighted variance
    weighted_variance = sum(weights .* euclidean_distances.^2) / sum(weights);
    disp(['Weighted Variance for ', power_variable, ': ', num2str(weighted_variance)]);

    % Compute the average power per patient and hemisphere
    [G, patient_id, Hemisphere] = findgroups(data.patient_id, data.Hemisphere);
    averaged_data = table(patient_id, Hemisphere, splitapply(@(x) mean(x, 'omitnan'), data.(power_variable), G), ...
        'VariableNames', {'patient_id', 'Hemisphere', [power_variable, '_per_p']});  % Correct column name

    % Join the averaged data with the original data
    joined = outerjoin(data, averaged_data, 'Keys', {'patient_id', 'Hemisphere'}, 'MergeKeys', true);
    joined.residual = joined.(power_variable) - joined.([power_variable, '_per_p']);  % Use dynamic reference

    % Initialize array for randomized variances
    random_weighted_variances = zeros(nperms, 1);

    for i = 1:nperms
        % Shuffle averaged power
        shuffled_avgs = averaged_data.([power_variable, '_per_p'])(randperm(height(averaged_data)));
        rand_avgs = averaged_data;
        rand_avgs.([power_variable, '_per_p']) = shuffled_avgs;

        % Merge with joined data
        rand_joined = outerjoin(joined, rand_avgs, 'Keys', {'patient_id', 'Hemisphere'}, 'MergeKeys', true);
        
        % Shuffle residuals within each patient
        [G, patient_id, Hemisphere] = findgroups(rand_joined.patient_id, rand_joined.Hemisphere);
        rand_joined.shuffled_residual = nan(height(rand_joined), 1);
        groupCell = splitapply(@(x) {x(randperm(numel(x)))}, rand_joined.residual, G);

        for g = 1:numel(groupCell)
            rand_joined.shuffled_residual(G == g) = groupCell{g};
        end

        % Compute randomized power
        rand_joined.rand_power_per_p = rand_joined.shuffled_residual + rand_joined.([power_variable, '_per_p_joined']);
        rand_weights = rand_joined.rand_power_per_p;

        % Compute randomized weighted averages
        rand_x_avg = sum(rand_joined.x_flip .* rand_weights) / sum(rand_weights);
        rand_y_avg = sum(rand_joined.y_flip .* rand_weights) / sum(rand_weights);
        rand_z_avg = sum(rand_joined.z_flip .* rand_weights) / sum(rand_weights);

        % Compute Euclidean distances
        rand_euclidean_distances = sqrt((rand_joined.x_flip - rand_x_avg).^2 + (rand_joined.y_flip - rand_y_avg).^2 + (rand_joined.z_flip - rand_z_avg).^2);

        % Compute weighted variance
        random_weighted_variances(i) = sum(rand_weights .* rand_euclidean_distances.^2) / sum(rand_weights);


    end

    % Compute p-value
    p_value = sum(random_weighted_variances <= weighted_variance) / (nperms + 1);
    disp(['P-value for ', power_variable, ': ', num2str(p_value)]);

    % Plot histogram for the current variable
    histogram(random_weighted_variances, floor(sqrt(nperms)), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'black');
    title(['Histogram of Inertial Variances for ', pow_label]);
    xlabel('Distance from Peak');
    ylabel('Frequency');
    hold on;

    % Annotate with weighted variance and p-value
    text(max(random_weighted_variances) * 0.99, max(ylim) * 0.8, ...
        sprintf('I = %.2f\np = %.3f', weighted_variance, p_value), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'red');

    % Add vertical line for weighted variance
    xline(weighted_variance, 'b', 'LineWidth', 2);
    hold off;

    % store stats
    stat = table();
    stat.Biomarker = string(Y_vars{idx});
    stat.wI = weighted_variance;
    stat.pwI = p_value;
    stat.rI = sqrt(weighted_variance);

    stats = [stats ;stat];
end



