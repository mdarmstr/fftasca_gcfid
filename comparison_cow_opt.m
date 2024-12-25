%% =====================================================================
%   COMPARISON SCRIPT: Optimize COW parameters for EACH jitter level
%
%   Outline:
%     1) Setup design, baseline "p_known"
%     2) For each Jitter Level:
%        (a) Generate ONE "representative" dataset
%        (b) Optimize (segment, slack) via optim_cow
%        (c) For each replicate:
%            - Generate data
%            - Raw analysis (time/frequency)
%            - COW alignment with best (segment, slack)
%            - COW analysis (time/frequency)
%     3) Convert p-values to Delta z-scores
%     4) Plot results
% =====================================================================

clear; close all; clc;

%% 1) BASIC SETTINGS
reps       = 4;        % For your design
vars       = 5;
levels     = {[1,2,3,4], [1,2,3]};
inner_reps = 5;        % how many replicates per jitter level you want
jitterLevels = 0:5:50; % example: step by 5, or 0:1:50 if you like

% Create your design matrix F for parglm
F = create_design(levels, 'Replicates', reps);

% We'll store p-values here:
p_time_all     = zeros(inner_reps, length(jitterLevels));
p_freq_all     = zeros(inner_reps, length(jitterLevels));
p_time_all_cow = zeros(inner_reps, length(jitterLevels));
p_freq_all_cow = zeros(inner_reps, length(jitterLevels));

%% 1A) Get "p_known" from some baseline data
X_base = zeros(size(F,1), vars);
for ii = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(ii));
    X_base(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + ...
                    repmat(randn(1, vars), length(idx), 1);
end

[~, parglmo] = parglm(X_base, F);
p_known      = parglmo.p(1);
z_known      = norminv(1 - p_known);

%% 2) LOOP OVER JITTER LEVELS
for jL = 1:length(jitterLevels)
    
    thisJitter = jitterLevels(jL);
    disp('--------------------------------------------------------');
    disp(['Jitter Level = ' num2str(thisJitter)]);
    disp('--------------------------------------------------------');
    
    % 2(a) Generate ONE "representative" dataset to optimize for this jitter
    %      We'll just take one replicate for that. 
    %      Then, we re-use the best (segment, slack) for all replicates.
    X_temp = zeros(size(F,1), vars);
    for ii = 1:length(levels{1})
        idx = find(F(:,1) == levels{1}(ii));
        X_temp(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + ...
                        repmat(randn(1, vars), length(idx), 1);
    end
    
    % Some random columns, scaling, permutation, etc.
    Xr   = randn(size(X_temp,1), 5);
    Xtmp = [X_temp, Xr];
    perm = randperm(size(Xtmp,2));
    Xtmp = Xtmp - min(Xtmp(:));
    Xtmp = Xtmp .* 5000;
    Xtmp = Xtmp(:, perm);

    % Create "representative" chromatographic data with "thisJitter"
    xAxis = linspace(0, 5000, 5000);
    allPeaks = zeros(size(Xtmp,1), length(xAxis), size(Xtmp,2));
    for iSamp = 1:size(Xtmp,1)
        for jPeak = 1:size(Xtmp,2)
            amp    = Xtmp(iSamp, jPeak);
            pos    = 500 + (jPeak-1)*thisJitter*3;
            jitt   = thisJitter * randn();
            posJit = pos + jitt;
            
            sigma  = 5;
            allPeaks(iSamp,:,jPeak) = amp * ...
                exp(-(xAxis - posJit).^2/(2*sigma^2));
        end
    end
    syntheticChromDataRep = sum(allPeaks, 3);
    syntheticChromDataRep = syntheticChromDataRep + ...
                            20.*randn(size(syntheticChromDataRep)) + 100;

    % 2(b) Optimize COW on that "representative" data
    %      We'll pick the first row as reference (or average).
    refRep = syntheticChromDataRep(1,:);
    % refRep = mean(syntheticChromDataRep,1);  % if you prefer average

    optim_space = [20 80 1 10];   % e.g. segment in [20..80], slack in [1..10]
    opt_options = [0 3 50 0.15];  % [plot=0, #gridMax=3, maxSteps=50, bandFrac=0.15]
    
    disp('Optimizing (segment, slack) for this jitter level...');
    [optim_pars, ~, ~] = optim_cow(syntheticChromDataRep, ...
                                   optim_space, ...
                                   opt_options, ...
                                   refRep);

    best_segment = optim_pars(1);
    best_slack   = optim_pars(2);
    disp(['  => Best segment=', num2str(best_segment), ...
          ', slack=', num2str(best_slack)]);

    % 2(c) Now do the usual replicate loop for this jitter level
    parfor r = 1:inner_reps
        
        % i) Generate random "design" data
        X = zeros(size(F,1), vars);
        for ii = 1:length(levels{1})
            idx = find(F(:,1) == levels{1}(ii));
            X(idx,:) = simuleMV(length(idx), vars, 'LevelCorr', 8) + ...
                       repmat(randn(1, vars), length(idx), 1);
        end
        
        Xr   = randn(size(X,1),5);
        X    = [X,Xr];
        perm = randperm(size(X,2));
        X    = X - min(X(:));
        X    = X.*5000;
        X    = X(:,perm);
        
        % Make chromatographic data with "thisJitter"
        nRows = size(X,1);
        nCols = size(X,2);
        allPeaks = zeros(nRows, length(xAxis), nCols);
        for iSamp = 1:nRows
            for jPeak = 1:nCols
                amp    = X(iSamp,jPeak);
                pos    = 500 + (jPeak-1)*thisJitter*3;
                jitt   = thisJitter*randn();  
                posJit = pos + jitt;
                sigma  = 5;
                allPeaks(iSamp,:,jPeak) = amp .* ...
                    exp(-(xAxis - posJit).^2/(2*sigma^2));
            end
        end
        
        syntheticChromData = sum(allPeaks,3);
        syntheticChromData = syntheticChromData + ...
                             20.*randn(size(syntheticChromData)) + 100;
        
        % --- RAW TIME
        [~, parglmot] = parglm(syntheticChromData, F, 'Preprocessing', 1);
        p_time_all(r,jL) = parglmot.p(1);

        % --- RAW FREQ
        fftChromData = fft(syntheticChromData, [], 2);
        halfN        = floor(size(fftChromData,2)/2);
        fftChromData = fftChromData(:,1:halfN+1);
        [~, parglmof] = parglm(fftChromData, F, 'Preprocessing', 1);
        p_freq_all(r,jL) = parglmof.p(1);

        % --- COW alignment w/ best_segment, best_slack
        ref = syntheticChromData(1,:);  % or average
        [~, syntheticChromDataCOW, ~] = cow(ref, ...
                                            syntheticChromData, ...
                                            best_segment, ...
                                            best_slack, ...
                                            [0 1 0 round(length(ref)*opt_options(4)) 0]);

        % --- COW TIME
        [~, parglmot_cow] = parglm(syntheticChromDataCOW, F, 'Preprocessing', 1);
        p_time_all_cow(r,jL) = parglmot_cow.p(1);

        % --- COW FREQ
        fftChromDataCOW = fft(syntheticChromDataCOW, [], 2);
        fftChromDataCOW = fftChromDataCOW(:,1:halfN+1);
        [~, parglmof_cow] = parglm(fftChromDataCOW, F, 'Preprocessing', 1);
        p_freq_all_cow(r,jL) = parglmof_cow.p(1);
        
    end % replicate loop
    
    disp(['   Done with all replicates at Jitter=' num2str(thisJitter)]);
end

%% 3) Convert p-values => z-scores (minus known baseline)
z_time_all     = norminv(1 - p_time_all)     - z_known;
z_freq_all     = norminv(1 - p_freq_all)     - z_known;
z_time_all_cow = norminv(1 - p_time_all_cow) - z_known;
z_freq_all_cow = norminv(1 - p_freq_all_cow) - z_known;

%% 4) Mean & std across replicates
z_time_mean     = mean(z_time_all,1);    z_time_std     = std(z_time_all,0,1);
z_freq_mean     = mean(z_freq_all,1);    z_freq_std     = std(z_freq_all,0,1);
z_time_mean_cow = mean(z_time_all_cow,1);z_time_std_cow = std(z_time_all_cow,0,1);
z_freq_mean_cow = mean(z_freq_all_cow,1);z_freq_std_cow = std(z_freq_all_cow,0,1);

%% 5) Plot results
figure; hold on;

% Raw freq (line+shaded)
freqColor = [0 0.35 0.9];
ln_raw_freq = plot(jitterLevels, z_freq_mean, 'Color', freqColor, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_freq_mean - z_freq_std, fliplr(z_freq_mean + z_freq_std)], ...
     freqColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Raw time (line+shaded)
timeColor = [0.9 0.2 0];
ln_raw_time = plot(jitterLevels, z_time_mean, 'Color', timeColor, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_time_mean - z_time_std, fliplr(z_time_mean + z_time_std)], ...
     timeColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% COW freq
freqColorCOW = [0 0.7 0.7];
ln_cow_freq = plot(jitterLevels, z_freq_mean_cow, '--', 'Color', freqColorCOW, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_freq_mean_cow - z_freq_std_cow, fliplr(z_freq_mean_cow + z_freq_std_cow)], ...
     freqColorCOW, 'FaceAlpha', 0.2, 'EdgeColor','none');

% COW time
timeColorCOW = [0.95 0.45 0];
ln_cow_time = plot(jitterLevels, z_time_mean_cow, '--', 'Color', timeColorCOW, 'LineWidth', 2);
fill([jitterLevels, fliplr(jitterLevels)], ...
     [z_time_mean_cow - z_time_std_cow, fliplr(z_time_mean_cow + z_time_std_cow)], ...
     timeColorCOW, 'FaceAlpha', 0.2, 'EdgeColor','none');

xlabel('Jitter Level');
ylabel('z_{obs} - z_{known}');
title('Comparison: Re-optimized COW for each Jitter Level');
hndl = [ln_raw_freq(1); ln_raw_time(1); ln_cow_time(1); ln_cow_freq(1)];
%legend({'Freq (raw)','Time (raw)','Freq (COW)','Time (COW)'}, 'Location','best');
legend(hndl,{'Freq (raw)','Time (raw)','Freq (COW)','Time (COW)'},'Location','best')
set(gca,'FontSize',14);

hold off;

% optionally export
exportgraphics(gcf, fullfile('figures','comparison_cow_reoptimized.pdf'), ...
               'ContentType','vector');

disp('All done!');
