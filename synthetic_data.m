reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};
ro = 5;


% Levels of jitter (you can modify the scale)
jitterLevels = 0:5:50;  % Five levels of jitter

p_time = zeros(1,length(jitterLevels));
p_freq = zeros(1,length(jitterLevels));

reps = 4;
vars = 5;
levels = {[1,2,3,4],[1,2,3]};

F = create_design(levels,'Replicates',reps);

X = zeros(size(F,1),vars);
for ii = 1:length(levels{1}),
  X(find(F(:,1) == levels{1}(ii)),:) = simuleMV(length(find(F(:,1) == levels{1}(ii))),vars,'LevelCorr',8) + repmat(randn(1,vars),length(find(F(:,1) == levels{1}(ii))),1);
end
 

Xr = randn(size(X,1),5);

X = [X,Xr];
perm = randperm(size(X,2));

X = X - min(min(X));

X = X.*5000;

X = X(:,perm);

% Define the synthetic data matrix
dataMatrix = X;

% Parameters for the Gaussian function
x = linspace(0, 5000, 5000);  % x-axis (chromatographic time scale)

% Preallocate a storage matrix for jittered data
jitteredPeaks = zeros(size(dataMatrix, 1), length(x), size(dataMatrix, 2));
for jitterLevel = 1:length(jitterLevels)
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
    
    % Sum the Gaussian peaks for each sample
    syntheticChromData = sum(jitteredPeaks, 3);
    syntheticChromData = syntheticChromData + 20.*randn(size(syntheticChromData,1),size(syntheticChromData,2)) + 100;
    
    % Plot the chromatographic data with jitter for this level
    figure;
    hold on;
    for i = 1:size(dataMatrix, 1)
        plot(syntheticChromData(i,1:100));
    end
    title(['Synethetic Data with Jitter Level ' num2str(jitterLevel)]);
    xlabel('Acquisitions');
    ylabel('Intensity');
    hold off;

    [~,parglmot] = parglm(syntheticChromData,F,'Model',{[1,2]},'Preprocessing',1);
    p_time(jitterLevel) = parglmot.p(1);

    fftChromData = fft(syntheticChromData,[],2);
    halfN = floor(size(fftChromData,2)/2);
    fftChromData = fftChromData(:,1:halfN+1);
    [~,parglmof] = parglm(fftChromData,F,'Model',{[1,2]},'Preprocessing',1);
    p_freq(jitterLevel) = parglmof.p(1);

    disp('working..')

end


plot(jitterLevels,p_freq);
hold on;
plot(jitterLevels,p_time);
hold off;