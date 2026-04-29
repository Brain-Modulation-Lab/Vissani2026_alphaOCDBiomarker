% Test inertial variance with simulated scenarios
clear; clc;

% ------------------------------------------------------------------------
% Paths and input files
% -------------------------------------------------------------------------

scriptPath = pwd;
FOLDERPATH_DEMOS = fullfile(scriptPath, "demos");

bss_file = fullfile(FOLDERPATH_DEMOS, ...
    'sub-OP1002_ses-postop01_percept-bss_hemi-left.mat');

coord_file = fullfile(FOLDERPATH_DEMOS, ...
    'sub-OP1002_ses-postop01_lead-dbs-mni-coords_hemi-left.csv');

% ------------------------------------------------------------------------
% Load data
% -------------------------------------------------------------------------

tmp = load(bss_file);
selDataRaw = tmp.BipolarDB;

selDataCoords = readtable(coord_file);

Nuclei = {'GPe','BNST','NucleusAccumbens'};

% ------------------------------------------------------------------------
% Consolidate repeated recordings
% -------------------------------------------------------------------------

cons_keys = [
    "patient_id", ...
    "session_id", ...
    "Hemisphere", ...
    "Sensing_ChannelbipP", ...
    "Sensing_ChannelbipN"
];

selData = consolidate_within_repetition(selDataRaw, cons_keys);

% ------------------------------------------------------------------------
% Select power metric
% -------------------------------------------------------------------------

power_metric = 'MaxPowerAlphaFlat';
power_metric_mono = [power_metric, 'Contr'];
power_metric_label = [power_metric_mono, ' [dB]'];

if contains(power_metric, 'Flat')
    fvec = 'FrequencyFlat';
else
    fvec = 'Frequency';
end

% ------------------------------------------------------------------------
% Compute pseudo-monopolar estimate and join coordinates
% -------------------------------------------------------------------------

MonoPatternEst = get_MONOPATTERN(selData, power_metric);

MonoPatternData = join( ...
    MonoPatternEst, ...
    selDataCoords, ...
    "Keys", ["patient_id", "Hemisphere", "Sensing_Channel"] ...
);

if ~isnumeric(MonoPatternData.Sensing_Channel)
    MonoPatternData.Sensing_Channel = double(string(MonoPatternData.Sensing_Channel));
end

% ------------------------------------------------------------------------
% Flip right hemisphere coordinates
% -------------------------------------------------------------------------

right_idx = strcmp(string(MonoPatternData.Hemisphere), 'Right');

xyz = MonoPatternData{:, {'x','y','z'}};
xyz_flip = xyz;

if any(right_idx)
    rh_flip = ea_flip_lr_nonlinear(xyz);
    xyz_flip(right_idx, :) = rh_flip(right_idx, :);
end

MonoPatternData.x_flip = xyz_flip(:, 1);
MonoPatternData.y_flip = xyz_flip(:, 2);
MonoPatternData.z_flip = xyz_flip(:, 3);

% ------------------------------------------------------------------------
% Real-data anatomical distance
% -------------------------------------------------------------------------

[~, fh_dist] = anatomicaldistance_ephysio_target( ...
    MonoPatternData, ...
    Nuclei, ...
    power_metric_mono, ...
    power_metric_label);

sgtitle(fh_dist, "Distance from Nucleus");

% ------------------------------------------------------------------------
% Simulate inertial-variance scenarios
% -------------------------------------------------------------------------
% ephysio_focality_inertia expects:
% patient_id, Hemisphere, x_flip, y_flip, z_flip, and one or more power vars.
% It computes weighted spatial inertia and compares it to a permutation null.
% See function structure in your uploaded code. 

% Settings
rng(1);

nperms     = 1000;
n_patients = 60;
hemispheres = ["Left", "Right"];
channels = 0:3;

power_vars = { ...
    'LowInertia_Focal', ...
    'HighInertia_Diffuse'};

power_labels = { ...
    'Low inertia: focal hotspot', ...
    'High inertia: diffuse power'};

% ------------------------------------------------------------------------
% Generate simulated contact-level dataset
% -------------------------------------------------------------------------

sim_data = table();
row_id = 0;

for p = 1:n_patients

    patient_id = sprintf("SIM%03d", p);
    patient_shift = randn(1, 3) .* [1.5 1.5 1.5];

    for h = 1:numel(hemispheres)

        hemi = hemispheres(h);

        if hemi == "Left"
            hemi_center = [-10 4 -3] + patient_shift;
        else
            hemi_center = [10 4 -3] + patient_shift;
        end

        for ch = channels

            row_id = row_id + 1;

            % Simulated contacts arranged approximately along a DBS lead
            sim_data.patient_id(row_id,1) = string(patient_id);
            sim_data.Hemisphere(row_id,1) = hemi;
            sim_data.Sensing_Channel(row_id,1) = ch;

            sim_data.x(row_id,1) = hemi_center(1) + randn * 0.8;
            sim_data.y(row_id,1) = hemi_center(2) + ch * 1.5 + randn * 0.6;
            sim_data.z(row_id,1) = hemi_center(3) + ch * 1.2 + randn * 0.6;
        end
    end
end

% Flip right hemisphere to common left-space
sim_data.x_flip = sim_data.x;
sim_data.y_flip = sim_data.y;
sim_data.z_flip = sim_data.z;

right_idx = sim_data.Hemisphere == "Right";
sim_data.x_flip(right_idx) = -sim_data.x(right_idx);

xyz = [sim_data.x_flip, sim_data.y_flip, sim_data.z_flip];
n = height(sim_data);

% ------------------------------------------------------------------------
% Define hotspots and distances
% -------------------------------------------------------------------------

global_hotspot = [-10 6 -2];

d_global = sqrt(sum((xyz - global_hotspot).^2, 2));



% ------------------------------------------------------------------------
% Scenario 1: LOW inertia, bilateral/global focal hotspot
% -------------------------------------------------------------------------

sigma_low = 1.2;

low_inertia = exp(-(d_global.^2) ./ (2 * sigma_low^2));
low_inertia = low_inertia + 0.04 * randn(n, 1);

% ------------------------------------------------------------------------
% Scenario 2: HIGH inertia, diffuse power across all contacts
% -------------------------------------------------------------------------

sigma_high = 55;

high_inertia = exp(-(d_global.^2) ./ (2 * sigma_high^2));
high_inertia = high_inertia + 0.50 * rand(n, 1);




% Make all weights positive
sim_data.LowInertia_Focal = make_positive(low_inertia);
sim_data.HighInertia_Diffuse = make_positive(high_inertia);



% ------------------------------------------------------------------------
% 3D visualization of simulated spatial power patterns
% -------------------------------------------------------------------------

nt  = numSubplots(numel(power_vars));
%
fh_3d = figure( ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [100 100  nt(1)*1500 nt(2)*200 ]);

t3 = tiledlayout(fh_3d, nt(1), nt(2), ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for ii = 1:numel(power_vars)

    ax = nexttile(t3, ii);
    hold(ax, 'on');

    vals = sim_data.(power_vars{ii});

    scatter3( ...
        ax, ...
        sim_data.x_flip, ...
        sim_data.y_flip, ...
        sim_data.z_flip, ...
        70, ...
        vals, ...
        'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 0.25, ...
        'MarkerFaceAlpha', 0.20);

    % Show weighted center
    w = vals;
    wc = [
        sum(sim_data.x_flip .* w) ./ sum(w), ...
        sum(sim_data.y_flip .* w) ./ sum(w), ...
        sum(sim_data.z_flip .* w) ./ sum(w)
    ];

    plot3(ax, wc(1), wc(2), wc(3), ...
        'kp', ...
        'MarkerSize', 18, ...
        'MarkerFaceColor', 'y', ...
        'LineWidth', 1.5);

    xlabel(ax, 'x flip');
    ylabel(ax, 'y flip');
    zlabel(ax, 'z flip');

    title(ax, power_labels{ii}, ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold');

    axis(ax, 'equal');
    grid(ax, 'on');
    view(ax, 35, 22);

    cb = colorbar(ax);
    cb.Label.String = 'Simulated power';

    set(ax, ...
        'FontName', 'Helvetica', ...
        'FontSize', 11, ...
        'TickDir', 'out', ...
        'LineWidth', 0.75, ...
        'Box', 'off');
end

sgtitle(t3, '3D simulated spatial power scenarios');


% ------------------------------------------------------------------------
% Run inertial variance analysis on simulated data only
% -------------------------------------------------------------------------

[stats_sim, fh_sim] = ephysio_focality_inertia( ...
    sim_data, ...
    power_vars, ...
    power_labels, ...
    nperms);

set(fh_sim, 'Color', 'w', 'Position', [100 100 1400 900]);
sgtitle(fh_sim, 'Inertial variance null distributions: simulated data only');

disp(stats_sim);