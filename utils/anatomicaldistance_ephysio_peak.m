function [stats, fh]=  anatomicaldistance_ephysio_peak(data, Y_vars,Y_labels)

% correlation with distance to 95th percentile peak and alpha power
fh = figure("position",[200 200 1200 900]);

stats = table();

% Loop through each variable in Y_vars
for idx = 1:length(Y_vars)
    nexttile
    power_variable = Y_vars{idx};  % Current variable
    power_label = Y_labels{idx};

    % Filter data based on the 95th percentile of the current variable
    quantile_val = quantile(data.(power_variable), 0.95);
    centroid_data = data(data.(power_variable) >= quantile_val, :);

    % Compute centroid coordinates for the current variable
    centroid = [mean(centroid_data.x_flip, 'omitnan'), mean(centroid_data.y_flip, 'omitnan'), mean(centroid_data.z_flip, 'omitnan')];

    % Add centroid back to the original dataset and compute distances
    centroid_x = centroid(1);
    centroid_y = centroid(2);
    centroid_z = centroid(3);

    distance_x = abs(data.x_flip - centroid_x);
    distance_y = abs(data.y_flip - centroid_y);
    distance_z = abs(data.z_flip - centroid_z);
    distance = sqrt(distance_x.^2 + distance_y.^2 + distance_z.^2);

    fprintf("Centroid of peak (95th prctile): %s \n", power_variable)
    fprintf("X = %1.2f, Y = %1.2f ; Z = %1.2f [mm] \n", centroid_x, centroid_y, centroid_z)

    % Perform correlation test between distance and the current variable
    %[R, P] = corr(data.distance, data.(power_variable), 'Type', 'Pearson', 'Rows', 'complete');
    [R, pR] = compare_stat_scatterplot(distance, data.(power_variable), 'corr_sp');  % Compute correlation

    % computr 95 ci
    [~,~,R_l,R_u] = corrcoef(tiedrank(distance), tiedrank(data.(power_variable))); R_l = R_l(1,2); R_u = R_u(1,2);
     




    % Plotting
    %scatter(distance, data.(power_variable), 'b'); % scatter plot
    scatter(distance, data.(power_variable), 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2)

    hold on;


    % Fit linear model for showing regression line
    mdl = fitlm(distance, data.(power_variable));

    % Generate predictions
    X_pred = linspace(min(distance), max(distance), 100)'; % Smooth X-axis points
    [Y_pred, Y_CI] = predict(mdl, X_pred); % Predicted values and confidence intervals

    % Plot regression line
    plot(X_pred, Y_pred, 'r', 'LineWidth', 2);

    % Confidence interval shading
    fill([X_pred; flipud(X_pred)], [Y_CI(:,1); flipud(Y_CI(:,2))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');



    % Add annotations
    title([power_label, ' vs. Distance from Peak']);
    xlabel('Distance from Peak');
    ylabel(power_variable);
    text(9, max(data.(power_variable)), sprintf('r = %.2f [%.2f %.2f]\np = %.3f', R, R_l, R_u, pR), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'r');

    % Reverse x-axis for visual similarity to ggplot
    set(gca, 'XDir', 'reverse');
    hold off;

    stat = table();
    stat.Biomarker = string(Y_vars{idx});
    stat.R = R;
    stat.R_l = R_l;
    stat.R_u = R_u;
    stat.pR = pR;

    stats = [stats ;stat];
end
% 
% % save current figure
% figname = strcat('sub-all','_device-all','_datatype-brainsensesurvey','_session-postop','_session-first','_montage-','bipolar','_anatomy-distance-centroid-95thpeak');
% saveFigures(gcf,char(fullfile(PATH_FIGURES,'first',figname)),{'.png','.svg','.fig'})
% 
