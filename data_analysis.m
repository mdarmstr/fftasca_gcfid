%% A frequency-domain analysis using ASCA
clear all
close all

%Add the MEDA toolbox
addpath('MEDA/')
addpath('data/')

%% Peak table labels

lbls = {'RI929','RI1015','Methyl-1,4-benzoquinone','RI1079','Ethyl-1,4-benzoquinone','RI1174', 'RI1180','RI1230','RI1257','RI1287','RI1393','RI1410','RI1456','RI1467','RI1473','1,6-C15-diene','1-C15-ene','RI1540','RI1566','RI1577','RI1587','1-C16-ene','RI1608','RI1669','1,8-C17-diene','1-C17-ene'}';


%% Colourmap

red = rgb2hsv([202,45,48]./255);
tel = rgb2hsv([0,140,186]./255);
blk = [52,53,54]./255;

col_mat_table = hsv2rgb(interpolateColors(red,tel,100));
col_mat_table = [[0,0,0]; col_mat_table];

%% Peak table analysis
Xp = readmatrix('X_peak.csv');
Yp = readmatrix('Y_peak.csv');

blk_idx = sum(Yp,2) == 0;

Xp(blk_idx,:) = [];
Yp(blk_idx,:) = [];

[Ttable1,parglmo] = parglm(Xp,Yp,'Model',{[1,2]},'Preprocessing',1);
save_parglm = parglmo;
%[Ttable2, ~] = parglm_cell2(Xp, Yp, [1,2;1,3;1,4;2,3;2,4;3,4],2);
[Ttable2,parglmo_cell] = parglm_cell2(Xp, Yp, [1,2],1);

plot_matrix(Xp,col_mat_table,'Observations','Variables','X_peak_block');
plot_matrix(parglmo.D,col_mat_table(2:end,:),'Observations','Design','Y_peak_block');
plot_matrix(parglmo.B,col_mat_table(2:end,:),' ','Coefficients','B_peak_block');
plot_matrix(parglmo.residuals,col_mat_table,'Observations','Variables','E_peak_block');

%[tblp,parglmop] = parglm(Xp,Yp,'Model',{[1,2],[1,4]},'Preprocessing',2,'Ts',2);

Ttable2.Source{2} = 'Time-1';
Ttable2.Source{3} = 'Treatment-2';
Ttable2.Source{4} = 'Sex-3';
Ttable2.Source{5} = 'Order-4';

ascao = asca(parglmo);

%Plot factor time
no_clrs = length(unique(Yp(:,1)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
plot_vecm(ascao.factors{1}.scoresV,'ObsClass',Yp(:,1),'Colormap',clr_map)
set(gca, 'XTickLabel', []);
ylabel('PC1','Color',blk);
title('Factor - Time','Color',blk);
box off;
legend('24h','72h')
exportgraphics(gcf, fullfile('figures', ['peak_time', '.eps']), 'BackgroundColor', 'none');
%Plot factor time
no_clrs = length(unique(Yp(:,2)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
plot_scatterm(ascao.factors{2}.scoresV(:,1:2),'ObsClass',Yp(:,2),'Colormap',clr_map,'BlurIndex',0)
xlabel('PC1 Scores','Color',blk)
ylabel('PC2 Scores','Color',blk);
title('Factor - Treatment','Color',blk);
box off;
legend('Healthy','Primed','Wounded')
exportgraphics(gcf, fullfile('figures', ['peak_treat', '.eps']), 'BackgroundColor', 'none');
close all;

table2latex(Ttable2,'tables/peak_table.tex','%.3g',-1)
table2latex(Ttable1,'tables/peak_table_zero.tex','%.3g',1)
disp(Ttable2)

%% Do the analysis in the frequency domain.

X = readmatrix('X_data.csv')';
Y = readmatrix('Y_data.csv');

blk_idx = sum(Y,2) == 0;

% Blank runs for the two time regimes
blkidx24 = 1:4;
blkidx72 = 53:56;

% Remove all blanks from the data
X(blk_idx,:) = [];
Y(blk_idx,:) = [];

% Split data into pre-processed and original frequency representation
Xo = fft(X,[],2);
Xp = Xo;

% Preprocessing the data
x72p = fft(X(blkidx72,:),[],2);
x24p = fft(X(blkidx24,:),[],2);

x72m = mean(x72p);
x72s = std(x72p);

x24m = mean(x24p);
x24s = std(x24p);

% Blank removal according to the levels in Y
Xp(Y(:,1)==1,:) = (Xp(Y(:,1) == 1,:) - x24m);%./x24s;
Xp(Y(:,1)==2,:) = (Xp(Y(:,1) == 2,:) - x72m);%./x72s;

Fs = 200;            % Sampling frequency                    
tp = 1/Fs;             % Sampling period       
L = 45001;             % Length of signal

% Split the Frequency spectrum in half
Xf = zeros(size(X,1), size(X,2));
halfN = floor(size(Xf,2)/2);
Xf1 = Xo(:,1:halfN+1);
Xf2 = Xp(:,1:halfN+1);

% [Xf1,mn1,std1] = preprocess2D(Xf1,'preprocessing',1);
% [Xf2,mn2,std2] = preprocess2D(Xf2,'preprocessing',2);

[t_freq_mn,parglmo] = parglm(Xf2,Y,'Model',{[1,2]},'Preprocessing',1);

line_map = hsv2rgb(interpolateColors(red,tel,size(Xf,1)));

hold on;
for ii = 1:size(Xf,1)
    plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Xo(ii,:)))+ii,'Color',[line_map(ii,:),0.6]);
end
hold off;

xlim([-5,5]);
xlabel('f(Hz)','Color',blk);
ylabel('Amplitude','Color',blk);
title('Frequency Domain Representation','Color',blk);
box off;
exportgraphics(gcf, fullfile('figures', ['freq_rep', '.eps']), 'BackgroundColor', 'white');
close all;

hold on;
for ii = 1:size(X,1)
    plot(X(ii,17400:18000),'Color',[line_map(ii,:),0.6]);
end
hold off;

xlabel('Acquisitions','Color',blk);
ylabel('Signal Amplitude','Color',blk);
title('Time Domain Representation','Color',blk);
box off;
exportgraphics(gcf, fullfile('figures', ['time_rep', '.eps']), 'BackgroundColor', 'white');
close all;

% Mean centering the data
X_treat = parglmo.factors{2}.matrix;
X_treat_flip = fliplr(Xf1(:,1:halfN));
X_treat = [X_treat, X_treat_flip];

X_treat = ifft(X_treat,[],2);

class1_indices = find(Y(:,2) == 1);
class2_indices = find(Y(:,2) == 2);
class3_indices = find(Y(:,2) == 3);

figure;
hold on;

% Plot a single representative line for each class for the legend
h1 = plot(real(X_treat(class1_indices(1), 13888:14512)), 'b', 'DisplayName', 'Healthy');
h2 = plot(real(X_treat(class2_indices(1), 13888:14512)), 'r', 'DisplayName', 'Primed');
h3 = plot(real(X_treat(class3_indices(1), 13888:14512)), 'k', 'DisplayName', 'Wounded');


% Plot class 1 in blue
for i = 1:length(class1_indices)
    plot(real(X_treat(class1_indices(i), 13888:14512)), 'b');
end

% Plot class 2 in red
for i = 1:length(class2_indices)
    plot(real(X_treat(class2_indices(i), 13888:14512)), 'r');
end

% Plot class 2 in black
for i = 1:length(class3_indices)
    plot(real(X_treat(class3_indices(i), 13888:14512)), 'k');
end

xlabel('Acquisitions');
ylabel('Amplitude, centered');
title('Treatment - reconstruction');
legend([h1,h2,h3]);

ax_h = get(gcf,'Children');

 for i = 1:length(ax_h)
    if strcmp(get(ax_h(i), 'type'), 'axes')
        set(ax_h(i), 'FontSize', 14);
        val = i;
    end
end

exportgraphics(gcf, fullfile('figures', ['time_treat', '.eps']), 'BackgroundColor', 'none');
close all;

t_freq_mn.Source{2} = 'Time-1';
t_freq_mn.Source{3} = 'Treatment-2';
t_freq_mn.Source{4} = 'Sex-3';
t_freq_mn.Source{5} = 'Order-4';

disp(t_freq_mn)
table2latex(t_freq_mn,'tables/freq_table_mn.tex','%.3g',2)

ascao = asca(parglmo);

%Plot factor time
%Scores
no_clrs = length(unique(Y(:,1)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
%plot_vecm(real(ascao.factors{1}.scoresV),'ObsClass',Y(:,1),'Colormap',clr_map)
plot_vecm(real(ascao.factors{1}.scoresV),'ObsClass',Y(:,1),'Colormap',clr_map)
%plot_scatterm([real(ascao.factors{3}.scoresV(:,1)),imag(ascao.factors{3}.scoresV(:,1))],'ObsClass',Y(:,3),'Colormap',clr_map)
set(gca, 'XTickLabel', []);
ylabel('PC1','Color',blk);
title('Factor Time, Frequency','Color',blk);
box off;
legend('24h','72h')
exportgraphics(gcf, fullfile('figures', ['freq_time_mn', '.eps']), 'BackgroundColor', 'none');
close all;

%Plot effect - loads
time_loads = ascao.factors{1}.loads;
time_loads_flip = (flipud(time_loads(1:halfN)));
time_loads_total = [time_loads_flip;time_loads];
N = length(time_loads_total);
time_loads_total(floor(N/2)+2:end) = conj(flipud(time_loads_total(2:ceil(N/2))));

ascao_peaks = asca(save_parglm);
pks_tbl = ascao_peaks.factors{1}.loads;

plot_vec(pks_tbl);
title('Peak table loadings, Factor Time')
xlabel('Position')
exportgraphics(gcf, fullfile('figures', ['peak_loads_time', '.eps']));

plot_complex_loadings(real(time_loads_total), 'Real frequency loadings, Time', Fs/L*(-L/2:L/2-1));
exportgraphics(gcf, fullfile('figures', ['freq_loads_time_real', '.eps']));

plot_complex_loadings(imag(time_loads_total), 'Imag Frequency loadings, Time', Fs/L*(-L/2:L/2-1));
exportgraphics(gcf, fullfile('figures', ['freq_loads_time_imag', '.eps']));


time_peaks = real(ifft(fftshift(time_loads_total)));
plot(time_peaks,'LineWidth',2)
xlabel('Acquisitions')
ylabel('Normalized amplitude')
title('Inverse loadings - Time')

ax_h = get(gcf,'Children');

for i = 1:length(ax_h)
    if strcmp(get(ax_h(i), 'type'), 'axes')
        set(ax_h(i), 'FontSize', 14);
        val = i;
    end
end

exportgraphics(gcf, fullfile('figures', ['time_peaks', '.eps']));

no_clrs = length(unique(Y(:,2)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
plot_scatterm(real(ascao.factors{2}.scoresV(:,1:2)),'ObsClass',Y(:,2),'Colormap',clr_map,'BlurIndex',0)
xlabel('PC1 Scores','Color',blk)
ylabel('PC2 Scores','Color',blk);
title('Factor Treatment, Autoscaled','Color',blk);
box off;
legend('Healthy','Primed','Wounded')

exportgraphics(gcf, fullfile('figures', ['freq_treat_mn', '.eps']), 'BackgroundColor', 'none');
close all;

% autoscaling
[t_freq_auto,parglmo] = parglm(Xf2,Y,'Model',{[1,2]},'Preprocessing',2,'Permutations',1000);

t_freq_auto.Source{2} = 'Time-1';
t_freq_auto.Source{3} = 'Treatment-2';
t_freq_auto.Source{4} = 'Sex-3';
t_freq_auto.Source{5} = 'Order-4';

disp(t_freq_auto)
table2latex(t_freq_auto,'tables/freq_table_auto.tex','%.3g',-1)

ascao = asca(parglmo);

%Plot factor time
no_clrs = length(unique(Y(:,1)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
plot_vecm(real(ascao.factors{1}.scoresV),'ObsClass',Y(:,1),'Colormap',clr_map)
set(gca, 'XTickLabel', []);
ylabel('PC1','Color',blk);
title('Factor Time, AS','Color',blk);
box off;
legend('24h','72h')
exportgraphics(gcf, fullfile('figures', ['freq_time_auto', '.eps']), 'BackgroundColor', 'none');

no_clrs = length(unique(Y(:,2)));
clr_map = hsv2rgb(interpolateColors(red,tel,no_clrs));
plot_scatterm(real(ascao.factors{2}.scoresV(:,1:2)),'ObsClass',Y(:,2),'Colormap',clr_map,'BlurIndex',0)
xlabel('PC1 Scores','Color',blk)
ylabel('PC2 Scores','Color',blk);
title('Factor Treatment, SNV','Color',blk);
box off;
legend('Healthy','Primed','Wounded')
exportgraphics(gcf, fullfile('figures', ['freq_treat_auto', '.eps']), 'BackgroundColor', 'none');
close all;

% plot_vec(real(ascao.factors{1}.scoresV),'ObsClass',Y(:,1))
% plot_scatter(real(ascao.factors{2}.scoresV),'ObsClass',Y(:,2))
% plot_vec(real(ascao.factors{3}.scoresV),'ObsClass',Y(:,3))
% plot_scatter(real(ascao.factors{4}.scoresV(:,1:2)),'ObsClass',Y(:,4))

function plot_complex_loadings(z, plotTitle, frequencyRange)

% Define the x-axis values (frequency range)
x = frequencyRange;

% Create the plot
figure;

% Plot the real component in black with semi-transparency
plot(x, z, 'Color', [0 0 0,0.5], 'LineWidth', 2); % 'k' is the color code for black with alpha

% Plot the imaginary component in purple with semi-transparency

% Add labels and legend
xlabel('Frequency');
ylabel('Value');
title(plotTitle);

xlim([-5,5])

ax_h = get(gcf,'Children');

for i = 1:length(ax_h)
    if strcmp(get(ax_h(i), 'type'), 'axes')
        set(ax_h(i), 'FontSize', 14);
        val = i;
    end
end

end

function interpolatedColors = interpolateColors(color1, color2, numPoints)
    % interpolateColors creates a linear interpolation between two colors.
    %
    % Inputs:
    %   color1 - 1x3 vector representing the first RGB color
    %   color2 - 1x3 vector representing the second RGB color
    %   numPoints - The number of points in the interpolation
    %
    % Output:
    %   interpolatedColors - numPoints x 3 matrix of interpolated RGB colors

    % Create a linspace vector for interpolation
    t = linspace(0, 1, numPoints); % numPoints from 0 to 1

    % Initialize the matrix to hold the interpolated colors
    interpolatedColors = zeros(length(t), 3);

    % Perform the interpolation
    for i = 1:length(t)
        interpolatedColors(i, :) = t(i) * color2 + (1 - t(i)) * color1;
    end

end

function plot_matrix(mat, col_map, is_y, is_x, filename)
    % Normalize the matrix for imagesc
    imagesc(mat ./ max(mat));
    
    % Set the custom colormap
    colormap(col_map);
    
    % Get the size of the matrix
    [numRows, numCols] = size(mat);

    % Calculate the appropriate figure size in pixels
    dpi = 15; % Dots per inch for screen display, adjust as necessary
    screenPixelsPerUnit = 1; % Each unit in the data corresponds to 1 pixel on the screen
    
    % Calculate figure dimensions in inches
    figWidthInches = numCols / dpi;
    figHeightInches = numRows / dpi;

    % Set the figure size to match the image size
    set(gcf, 'Units', 'inches', 'Position', [1, 1, figWidthInches, figHeightInches]);

    % Add labels with larger, bold text
    xlabel(is_x, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
    ylabel(is_y, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

    % Remove x and y axis tick labels
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);

    % Position the x label at the top
    set(gca, 'XAxisLocation', 'top');

    % Remove the border
    set(gca, 'box', 'off');
    
    % Set the aspect ratio based on matrix dimensions
    %pbaspect([numCols numRows 1]);
    
    % Ensure the 'figures' directory exists
    if ~exist('figures', 'dir')
        mkdir('figures');
    end

    % Save the figure
    exportgraphics(gcf, fullfile('figures', [filename, '.eps']), 'BackgroundColor', 'none');

    close all;
end
