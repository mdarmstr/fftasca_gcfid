close all;

addpath('MEDA/')
addpath('WarpingTB/')

reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};
ro = 5;

% Levels of jitter (you can modify the scale)
jitterLevels = 0:5:50;  % Five levels of jitter

p_time = zeros(1,length(jitterLevels));
p_freq = zeros(1,length(jitterLevels));

% For storing the post-COW p-values
p_time_cow = zeros(1,length(jitterLevels));
p_freq_cow = zeros(1,length(jitterLevels));

reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};

F = create_design(levels,'Replicates',reps);

X = zeros(size(F,1),vars);
for ii = 1:length(levels{1})
    idx = find(F(:,1) == levels{1}(ii));
    % Insert your code that simulates correlated data
    X(idx,:) = simuleMV(length(idx),vars,'LevelCorr',8) + ...
               repmat(randn(1,vars),length(idx),1);
end

Xr = randn(size(X,1),5);
X = [X,Xr];
perm = randperm(size(X,2));
X = X - min(min(X));
X = X.*5000;
X = X(:,perm);

% Synthetic data matrix
dataMatrix = X;

% Parameters for the Gaussian function
x = linspace(0, 5000, 5000);  % x-axis (chromatographic time scale)

% COW parameters (adjust to your case)
Seg   = 50;                % Segment length (or boundaries). Tune as needed.
Slack = 5;                 % Slack. Tune as needed.
Opts  = [0 1 0 0 0];       % Options for cow [Plot, corrPower, fixSegLen, bandConstraint, saveDiagnostics]
%   Example: [0 1 0 0 0] means: no plotting, 1st power, no forced equal seg-len, no band constraint, no table stored

% Preallocate a storage matrix for jittered data
jitteredPeaks = zeros(size(dataMatrix, 1), length(x), size(dataMatrix, 2));

for jitterLevel = 1:length(jitterLevels)
    
    %---------------------------------------------------
    % 1) Generate new syntheticChromData for this jitter level
    %---------------------------------------------------
    for i = 1:size(dataMatrix, 1)  % Loop through each sample (row)
        for j = 1:size(dataMatrix, 2)  % Loop through each peak (column)
            amplitude = dataMatrix(i, j);
            peakPosition = 500 + (j-1)*jitterLevels(jitterLevel)*3;  % Assume fixed intervals between peaks
            
            % Introduce random jitter to the peak position
            jitter = jitterLevels(jitterLevel) * randn();
            jitteredPosition = peakPosition + jitter;
            
            % Create the Gaussian peak
            sigma = 5;  % Standard deviation of the Gaussian peak
            gaussianPeak = amplitude * exp(-(x - jitteredPosition).^2 / (2*sigma^2));
            
            % Store the jittered Gaussian peak
            jitteredPeaks(i,:,j) = gaussianPeak;
        end
    end
    
    syntheticChromData = sum(jitteredPeaks, 3);
    syntheticChromData = syntheticChromData + 20.*randn(size(syntheticChromData,1),size(syntheticChromData,2)) + 100;

    %---------------------------------------------------
    % 2) Plot the first few time profiles to visualize jitter
    %---------------------------------------------------
    figure;
    hold on;
    for i = 1:3
        plot(syntheticChromData(i,1:2000)');
    end
    
    title(['Synethetic Data with Jitter Level ' num2str(jitterLevel)]);
    xlabel('Acquisitions');
    ylabel('Intensity');
    hold off;

    ax_h = get(gcf,'Children');
    for i = 1:length(ax_h)
        if strcmp(get(ax_h(i), 'type'), 'axes')
            set(ax_h(i), 'FontSize', 14);
        end
    end

    exportgraphics(gcf, ...
        strcat("figures/jitter_illustration",sprintf('%d',jitterLevel),".pdf"), ...
        'ContentType','vector');
    
    %---------------------------------------------------
    % 3) GLM on raw data (no warping) 
    %---------------------------------------------------
    % a) Time domain
    [~,parglmot] = parglm(syntheticChromData,F,'Model',{[1,2]},'Preprocessing',1);
    p_time(jitterLevel) = parglmot.p(1);

    % b) Frequency domain
    fftChromData = fft(syntheticChromData,[],2);
    halfN        = floor(size(fftChromData,2)/2);
    fftChromData = fftChromData(:,1:halfN+1);
    
    [~,parglmof] = parglm(fftChromData,F,'Model',{[1,2]},'Preprocessing',1);
    p_freq(jitterLevel) = parglmof.p(1);
    
    %---------------------------------------------------
    % 4) COW Alignment of the time-domain data
    %---------------------------------------------------
    %  Choose the first sample as the reference T. 
    %  (Alternatively, use the mean or median across rows as your reference.)
    T = syntheticChromData(1,:);  
    
    %  Warp the entire matrix in one call:
    %  In cow.m, T must be [1 x pT], X must be [mP x nP].
    %  Here, mP = number of samples, nP = number of points per sample.
    %  So we pass T as above, and X as syntheticChromData.
    [~, syntheticChromDataCOW, ~] = cow(T, syntheticChromData, Seg, Slack, Opts);
    %  syntheticChromDataCOW will be [mP x length(T)] aligned to T.

    %---------------------------------------------------
    % 5) GLM on COW-aligned time series
    %---------------------------------------------------
    [~, parglmot_cow] = parglm(syntheticChromDataCOW,F,'Model',{[1,2]},'Preprocessing',1);
    p_time_cow(jitterLevel) = parglmot_cow.p(1);

    %---------------------------------------------------
    % 6) GLM on COW-aligned frequency domain
    %---------------------------------------------------
    fftChromDataCOW = fft(syntheticChromDataCOW,[],2);
    fftChromDataCOW = fftChromDataCOW(:,1:halfN+1);  % same half as before
    
    [~, parglmof_cow] = parglm(fftChromDataCOW,F,'Model',{[1,2]},'Preprocessing',1);
    p_freq_cow(jitterLevel) = parglmof_cow.p(1);
    
    disp('working..')

end

%-------------------------------------------------------
% 7) Compare p-values with/without warping
%-------------------------------------------------------
figure;
plot(jitterLevels,p_time,'-ob','LineWidth',1.5); 
hold on;
plot(jitterLevels,p_time_cow,'--ob','LineWidth',1.5); 
plot(jitterLevels,p_freq,'-or','LineWidth',1.5);
plot(jitterLevels,p_freq_cow,'--or','LineWidth',1.5);
xlabel('Jitter Level');
ylabel('p-value');
title('Effect of Jitter on p-values (Time vs. Frequency, with/without COW)');
legend({'Time (raw)','Time (COW)','Freq (raw)','Freq (COW)'},'Location','best');
set(gca,'FontSize',14);
hold off;
