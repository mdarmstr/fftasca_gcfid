reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};
inner_reps = 10;  % 5 replicates for each jitter level

% Levels of jitter
jitterLevels = 0:1:50;  % Five levels of jitter

% Storage for p-values across replicates and jitter levels
p_time_all = zeros(inner_reps, length(jitterLevels));
p_freq_all = zeros(inner_reps, length(jitterLevels));

F = create_design(levels,'Replicates',reps);

X = zeros(size(F,1),vars);
for ii = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(ii));
    X(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + repmat(randn(1,vars), length(idx), 1);
end

[~,parglmo] = parglm(X,F);

p_known = parglmo.p(1);

for jitterLevel = 1:length(jitterLevels)
    for r = 1:inner_reps  % 5 replicates for each jitter level
        X = zeros(size(F,1),vars);
        for ii = 1:length(levels{1})
            idx = find(F(:,1) == levels{1}(ii));
            X(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + repmat(randn(1,vars), length(idx), 1);
        end

        Xr = randn(size(X,1),5);
        X = [X,Xr];
        perm = randperm(size(X,2));
        X = X - min(min(X));
        X = X.*5000;
        X = X(:,perm);

        % Define the synthetic data matrix
        dataMatrix = X;
        x = linspace(0, 5000, 5000);  % x-axis (chromatographic time scale)

        % Preallocate a storage matrix for jittered data
        jitteredPeaks = zeros(size(dataMatrix, 1), length(x), size(dataMatrix, 2));

        for i = 1:size(dataMatrix, 1)  % Loop through each sample (row)
            for j = 1:size(dataMatrix, 2)  % Loop through each peak (column)
                amplitude = dataMatrix(i, j);
                peakPosition = 500 + (j-1) * jitterLevels(jitterLevel) * 3;  % Fixed intervals between peaks
                jitter = jitterLevels(jitterLevel) * randn();  % Random jitter for each peak
                jitteredPosition = peakPosition + jitter;

                % Create the Gaussian peak
                sigma = 5;  % Standard deviation of the Gaussian peak
                gaussianPeak = amplitude * exp(-(x - jitteredPosition).^2 / (2 * sigma^2));

                % Store the jittered Gaussian peak
                jitteredPeaks(i,:,j) = gaussianPeak;
            end
        end

        % Sum the Gaussian peaks for each sample
        syntheticChromData = sum(jitteredPeaks, 3);
        syntheticChromData = syntheticChromData + 20 .* randn(size(syntheticChromData)) + 100;

        % Compute p-values for time domain
        [~,parglmot] = parglm(syntheticChromData, F, 'Preprocessing', 1);
        p_time_all(r, jitterLevel) = parglmot.p(1);

        % Compute p-values for frequency domain
        fftChromData = fft(syntheticChromData, [], 2);
        halfN = floor(size(fftChromData, 2) / 2);
        fftChromData = fftChromData(:, 1:halfN+1);
        [~,parglmof] = parglm(fftChromData, F, 'Preprocessing', 1);
        p_freq_all(r, jitterLevel) = parglmof.p(1);
    end

    disp(['Working... Jitter Level: ' num2str(jitterLevels(jitterLevel))]);
end

% Compute mean and standard deviation for p-values across replicates
% p_time_mean = mean(p_time_all, 1);
% p_time_std = std(p_time_all, 0, 1);
% p_freq_mean = mean(p_freq_all, 1);
% p_freq_std = std(p_freq_all, 0, 1);

z_freq_all = norminv(1-p_freq_all);
z_time_all = norminv(1-p_time_all);
z_known    = norminv(1-p_known);

z_freq_all = z_freq_all - z_known;
z_time_all = z_time_all - z_known;

z_time_mean = mean(z_time_all,1);
z_time_std  = std(z_time_all,0,1);
z_freq_mean = mean(z_freq_all,1);
z_freq_std  = std(z_freq_all,0,1);


% Visualization
figure;
hold on;

% Line plot for mean p-values
plot(jitterLevels, z_freq_mean, 'b', 'LineWidth', 2);
plot(jitterLevels, z_time_mean, 'r', 'LineWidth', 2);

% Area plot for uncertainty (standard deviation)
fill([jitterLevels, fliplr(jitterLevels)], ...
    [z_freq_mean - z_freq_std, fliplr(z_freq_mean + z_freq_std)], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([jitterLevels, fliplr(jitterLevels)], ...
    [z_time_mean - z_time_std, fliplr(z_time_mean + z_time_std)], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xlabel('Jitter Levels');
ylabel('z_{obs} - z_{known}');
title('Comparison of analyses on synthetic data');
legend('\Delta z-values (frequency)', '\Delta z-values (time)', 'Uncertainty (frequency)', 'Uncertainty (time)', 'Location', 'best');
hold off;

ax_h = get(gcf,'Children');

 for i = 1:length(ax_h)
    if strcmp(get(ax_h(i), 'type'), 'axes')
        set(ax_h(i), 'FontSize', 14);
        val = i;
    end
 end

 exportgraphics(gcf, fullfile('figures', ['comparison', '.pdf']));

