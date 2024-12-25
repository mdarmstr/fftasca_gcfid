%% ----------------------------------------------------------------------
%  Main script with COW integration
%% ----------------------------------------------------------------------

clear; close all; clc;

addpath('MEDA/');
addpath('WarpingTB/')

reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};
inner_reps = 10;  % 10 replicates per jitter level

% Levels of jitter
jitterLevels = 0:1:50;  

% Storage for p-values across replicates and jitter levels (no COW)
p_time_all = zeros(inner_reps, length(jitterLevels));
p_freq_all = zeros(inner_reps, length(jitterLevels));

% Storage for p-values across replicates and jitter levels (with COW)
p_time_all_cow = zeros(inner_reps, length(jitterLevels));
p_freq_all_cow = zeros(inner_reps, length(jitterLevels));

% Setup design F
F = create_design(levels, 'Replicates', reps);

% Generate a baseline dataset just to get the "known" p-value
X_base = zeros(size(F,1), vars);
for ii = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(ii));
    X_base(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + ...
                    repmat(randn(1, vars), length(idx), 1);
end
[~, parglmo] = parglm(X_base, F);
p_known = parglmo.p(1);

%%
% COW parameters (tune as needed)
Seg   = 50;             % Segment length (or array of boundaries)
Slack = 5;              % Slack
Opts  = [0 1 0 0 0];    % COW options: [Plot=0, CorrPower=1, fixSegLen=0, bandConstraint=0, tableDiagnostics=0]

%%
% Main loop over jitter levels
parfor jitterLevel = 1:length(jitterLevels)
    
    for r = 1:inner_reps  
        
        % ----------------------------------------------------------
        % 1) Generate random data matrix X for the design
        % ----------------------------------------------------------
        X = zeros(size(F,1), vars);
        for ii = 1:length(levels{1})
            idx = find(F(:,1) == levels{1}(ii));
            X(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + ...
                       repmat(randn(1, vars), length(idx), 1);
        end
        
        % Some extra random columns, random permutation, scaling...
        Xr = randn(size(X,1),5);
        X  = [X, Xr];
        perm = randperm(size(X,2));
        X = X - min(X(:));
        X = X .* 5000;
        X = X(:,perm);
        
        % ----------------------------------------------------------
        % 2) Create synthetic chromatographic data with jitter
        % ----------------------------------------------------------
        dataMatrix = X;
        x = linspace(0, 5000, 5000);  
        
        jitteredPeaks = zeros(size(dataMatrix, 1), length(x), size(dataMatrix, 2));
        
        for iSamp = 1:size(dataMatrix, 1)
            for jPeak = 1:size(dataMatrix, 2)
                amplitude = dataMatrix(iSamp, jPeak);
                peakPosition = 500 + (jPeak-1) * jitterLevels(jitterLevel) * 3;
                jitter = jitterLevels(jitterLevel) * randn();  
                jitteredPosition = peakPosition + jitter;
                
                sigma = 5;  
                gaussianPeak = amplitude * exp(-(x - jitteredPosition).^2 / (2*sigma^2));
                
                jitteredPeaks(iSamp,:,jPeak) = gaussianPeak;
            end
        end
        
        syntheticChromData = sum(jitteredPeaks, 3);
        syntheticChromData = syntheticChromData + 20 .* randn(size(syntheticChromData)) + 100;
        
        % ----------------------------------------------------------
        % 3) GLM on time-domain (no COW)
        % ----------------------------------------------------------
        [~, parglmot] = parglm(syntheticChromData, F, 'Preprocessing', 1);
        p_time_all(r, jitterLevel) = parglmot.p(1);

        % ----------------------------------------------------------
        % 4) GLM on frequency-domain (no COW)
        % ----------------------------------------------------------
        fftChromData = fft(syntheticChromData, [], 2);
        halfN = floor(size(fftChromData, 2) / 2);
        fftChromData = fftChromData(:, 1:halfN+1);
        
        [~, parglmof] = parglm(fftChromData, F, 'Preprocessing', 1);
        p_freq_all(r, jitterLevel) = parglmof.p(1);

        % ----------------------------------------------------------
        % 5) COW Alignment (time-domain)
        % ----------------------------------------------------------
        % Choose sample #1 as the reference, or an average:
        T = syntheticChromData(1,:);  
        % T = mean(syntheticChromData,1);  % <--- Alternatively, use an average if you prefer

        % Warping each row to T
        %  cow() => [Warping, XWarped, Diagnos] = cow(T, X, Seg, Slack, Options)
        [~, syntheticChromDataCOW, ~] = cow(T, syntheticChromData, Seg, Slack, Opts);

        % ----------------------------------------------------------
        % 6) GLM on COW-aligned time-domain
        % ----------------------------------------------------------
        [~, parglmot_cow] = parglm(syntheticChromDataCOW, F, 'Preprocessing', 1);
        p_time_all_cow(r, jitterLevel) = parglmot_cow.p(1);

        % ----------------------------------------------------------
        % 7) GLM on COW-aligned frequency-domain
        % ----------------------------------------------------------
        fftChromDataCOW = fft(syntheticChromDataCOW, [], 2);
        fftChromDataCOW = fftChromDataCOW(:, 1:halfN+1);
        
        [~, parglmof_cow] = parglm(fftChromDataCOW, F, 'Preprocessing', 1);
        p_freq_all_cow(r, jitterLevel) = parglmof_cow.p(1);
        
    end
    
    disp(['Working... Jitter Level: ' num2str(jitterLevels(jitterLevel))]);
end

%%
% -----------------------------------------------------------------
% 8) Convert p-values to z-scores relative to "p_known" reference
% -----------------------------------------------------------------
z_known = norminv(1 - p_known);

% Raw data
z_time_all = norminv(1 - p_time_all) - z_known;
z_freq_all = norminv(1 - p_freq_all) - z_known;

% COW-aligned data
z_time_all_cow = norminv(1 - p_time_all_cow) - z_known;
z_freq_all_cow = norminv(1 - p_freq_all_cow) - z_known;

% Mean and Std across replicates
z_time_mean = mean(z_time_all, 1);
z_time_std  = std(z_time_all, 0, 1);
z_freq_mean = mean(z_freq_all, 1);
z_freq_std  = std(z_freq_all, 0, 1);

z_time_mean_cow = mean(z_time_all_cow, 1);
z_time_std_cow  = std(z_time_all_cow, 0, 1);
z_freq_mean_cow = mean(z_freq_all_cow, 1);
z_freq_std_cow  = std(z_freq_all_cow, 0, 1);

%%
% -----------------------------------------------------------------
% 9) Visualization of results (z-scores +/- std)
% -----------------------------------------------------------------
figure; hold on;

% --- (A) Raw frequency results
freqColor = [0 0.35 0.9];  % a shade of blue
timeColor = [0.9 0.2 0];

p1 = plot(jitterLevels, z_freq_mean, '-', 'Color', freqColor, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_freq_mean - z_freq_std, fliplr(z_freq_mean + z_freq_std)], ...
     freqColor, 'FaceAlpha', 0.2, 'EdgeColor','none');

% --- (B) Raw time results
p2 = plot(jitterLevels, z_time_mean, '-', 'Color', timeColor, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_time_mean - z_time_std, fliplr(z_time_mean + z_time_std)], ...
     timeColor, 'FaceAlpha', 0.2, 'EdgeColor','none');

% --- (C) COW-aligned frequency results
freqColorCOW = [0 0.7 0.7];  % teal-like
p3 = plot(jitterLevels, z_freq_mean_cow, '--', 'Color', freqColorCOW, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_freq_mean_cow - z_freq_std_cow, fliplr(z_freq_mean_cow + z_freq_std_cow)], ...
     freqColorCOW, 'FaceAlpha', 0.2, 'EdgeColor','none');

% --- (D) COW-aligned time results
timeColorCOW = [0.95 0.45 0];
p4 = plot(jitterLevels, z_time_mean_cow, '--', 'Color', timeColorCOW, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_time_mean_cow - z_time_std_cow, fliplr(z_time_mean_cow + z_time_std_cow)], ...
     timeColorCOW, 'FaceAlpha', 0.2, 'EdgeColor','none');

xlabel('Jitter Levels');
ylabel('z_{obs} - z_{known}');
title('Comparison of analyses (Time/Frequency) with/without COW');
legend([p1, p2, p3, p4], ...
       {'Freq (raw)', 'Time (raw)', 'Freq (COW)', 'Time (COW)'}, ...
       'Location','best');
set(gca, 'FontSize', 14);
hold off;

% Optional: export figure
exportgraphics(gcf, fullfile('figures', 'comparison_cow.pdf'), 'ContentType','vector');
