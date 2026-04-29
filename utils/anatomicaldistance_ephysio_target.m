function [stats, fh]=  anatomicaldistance_ephysio_target(Coords_atlas, Nuclei, Y_vars,Y_labels)


if ischar(Y_vars); Y_vars = {Y_vars}; end
if ischar(Y_labels); Y_labels = {Y_labels}; end


stats = table();
%Coords_atlas = Coords.(atlas);
fh = figure("Position",[200 100 300*numel(Y_vars) 300*numel(Nuclei)]);
tiledlayout(numel(Nuclei),numel(Y_vars))
%sgtitle(sprintf("Distance from Nucleus [%s]", atlas))

for nucl = 1 : numel(Nuclei)
    DataX = Coords_atlas.(Nuclei{nucl});  % X-axis data

    % Define Y-axis variables to be analyzed

    for i = 1:length(Y_vars)
        nexttile
        DataY = Coords_atlas.(Y_vars{i});  % Y-axis data
        [R, pR] = compare_stat_scatterplot(DataX, DataY, 'corr_sp');  % Compute correlation

        % Scatter plot
        scatter(DataX, DataY, 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 2)
        hold on;

        if numel(DataY) > 10
            % Fit linear model
            mdl = fitlm(DataX, DataY);

            % Generate predictions
            X_pred = linspace(min(DataX), max(DataX), 100)'; % Smooth X-axis points
            [Y_pred, Y_CI] = predict(mdl, X_pred); % Predicted values and confidence intervals

            % Plot regression line
            plot(X_pred, Y_pred, 'r', 'LineWidth', 2);

            % Confidence interval shading
            fill([X_pred; flipud(X_pred)], [Y_CI(:,1); flipud(Y_CI(:,2))], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');


            % computr 95 ci
            [~,~,R_l,R_u] = corrcoef(tiedrank(DataX), tiedrank(DataY)); R_l = R_l(1,2); R_u = R_u(1,2);
            stitle(sprintf("R = %1.2f [%.2f %.2f] [p = %1.4f]", R, R_l, R_u, pR))

        end


        hold off;
        xlabel(sprintf("Distance from %s", Nuclei{nucl}))
        ylabel(Y_labels{i})

        stat = table();

        if numel(DataY) > 10
            stat.Target = string(Nuclei{nucl});
            stat.Biomarker = string(Y_vars{i});
            stat.R = R;
            stat.R_l = R_l;
            stat.R_u = R_u;
            stat.pR = pR;
        end
        stats = [stats ;stat];
    end
end


