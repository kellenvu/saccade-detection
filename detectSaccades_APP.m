function detectSaccades_APP(params)
%{
Code by Kellen Vu.

Requirements
- Put .smr file in this directory to graph magnet data
- Put cali.mat file in this directory to scale magnet data
%}

%% Options
isUsingSaccades = params.isUsingSaccades;

isUsingCali = params.isUsingCali;
customScaleCh1 = params.customScaleCh1;
customScaleCh2 = params.customScaleCh2;

isRemovingSegment = params.isRemovingSegment;
removeTime = params.removeTime;

%% Load files

% Check mat files
matFiles = dir('*.mat*');
for i = 1:length(matFiles)
    if contains(matFiles(i).name, 'cali.mat')
        caliFile = matFiles(i);
    elseif contains(matFiles(i).name, 'saccades.mat')
        saccadeFile = matFiles(i);
    end
end
if isUsingCali && ~exist('caliFile', 'var')
    error('Make sure this .m file is in the same directory as the cali.mat file');
end

% Check for .smr file
if isempty(dir('*.smr'))
    error('Make sure this .m file is in the same directory as the .smr file');
end

% Load mat files
if isUsingCali
    load(caliFile.name, 'scaleCh1', 'scaleCh2');
end
if isUsingSaccades
    load(saccadeFile.name);
end

% Load data from .smr file
smrFile = dir('*.smr');
chanlist = readSpikeFile(smrFile.name, []);
chanindsAll = [chanlist.number];
chanlabels = {'hhpos', 'htpos', 'hepos1', 'hepos2', 'hepos', 'vepos', 'hhvel', 'htvel', 'htvel', 'TTL3', 'TTL4'};
chaninds = find(arrayfun(@(x) any(strcmp(x.title,chanlabels)), chanlist));
data = importSpike(smrFile.name, chanindsAll(chaninds));
if isUsingCali
    if scaleCh1 ~= 0
        magnetPos = scaleCh1 * datchandata(data, 'hepos1');
    else
        magnetPos = scaleCh2  * datchandata(data, 'hepos2');
    end
else
    magnetPos = customScaleCh1 * datchandata(data, 'hepos1') + customScaleCh2 * datchandata(data, 'hepos2');
end

%% Preprocess Magnet Data

magnetPos = double(magnetPos); % Convert to double
magnetPos = detrend(magnetPos); % Detrend
magnetPos = removeTransients(magnetPos, 1); % Remove transients

% Lowpass
N = 4;
fc = 40;
nq = 1000 / 2; % Half of sample rate
[bb,aa] = butter(N, fc / nq, 'low');
magnetPos = filtfilt(bb, aa, magnetPos);

%% Plotting

% Prep
filename = erase(smrFile.name, '.smr');
fig1 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
tiledlayout(2, 1);
h = sgtitle(filename);
h.Interpreter = 'none';
colormap = brewermap(4, 'Set1');

% Timescale
% magnetPos = magnetPos(0 * 1000 + 1:100 * 1000); % Crop timescale
magnetTime = (0:0.001:length(magnetPos) / 1000 - 0.001)';

% Magnet pos
ax1 = nexttile;
hold on;
plot(magnetTime, magnetPos, 'Color', [0.7, 0.7, 0.7]);

% Magnet vel
ax2 = nexttile;
hold on;
magnetVel = movingslope(magnetPos, 30) * 1000;
plot(magnetTime, magnetVel, 'Color', [0.7, 0.7, 0.7]);

% Detect saccades
if ~isUsingSaccades
    % Saccades
    saccadeCenters = detectSaccades(magnetTime, magnetPos, magnetVel, 1000, params.velThresFactor, true);
    saccades = getStartStop(saccadeCenters);
    
    if isRemovingSegment
        removeIndices = saccades(:, 1) > removeTime(1) * 1000 & saccades(:, 2) < removeTime(2) * 1000;
        saccades(removeIndices, :) = [];
    end
    
    % Get min PV
    peakVel = zeros(size(saccades, 1), 1);
    for i = 1:size(saccades, 1)
        start = saccades(i, 1);
        stop = saccades(i, 2);

        % Peak velocity
        absVel = abs(magnetVel);
        A = absVel(start:stop);
        peakVel(i) = max(A);
        minPV = min(peakVel);
    end
    
    saccades = extendSaccadesToThres(saccades, magnetPos, minPV);
    [saccades, artifacts] = removeArtifacts(saccades, magnetPos, params.maxArtifactDuration, params.neighborThres);
end

grid on;
xlabel('Time (s)');
ylabel('Velocity (deg / s)');
legend('Vel', 'Thres', 'Thres');
h = title("Velocity");
h.Interpreter = 'none';

% Plot
axes(ax1);
saccadesInContext = putInContext(saccades, magnetPos, true);
plot(magnetTime, saccadesInContext, 'Color', colormap(2, :), 'LineWidth', 1);
artifactsInContext = putInContext(artifacts, magnetPos, true);
plot(magnetTime, artifactsInContext, 'Color', colormap(1, :), 'LineWidth', 1);

grid on;
xlabel('Time (s)');
ylabel('Position (deg)');
h = title('Position');
h.Interpreter = 'none';
linkaxes([ax1, ax2], 'x');
legend('Pos', 'Saccades', 'Artifacts')

%% Analysis

% Extract features
duration = zeros(size(saccades, 1), 1);
peakVel = zeros(size(saccades, 1), 1);
amplitude = zeros(size(saccades, 1), 1);

for i = 1:size(saccades, 1)
    start = saccades(i, 1);
    stop = saccades(i, 2);
    
    % Peak velocity
    absVel = abs(magnetVel);
    A = absVel(start:stop);
    peakVel(i) = max(A);
    
    % Amplitude
    A = magnetPos(start:stop);
    amplitude(i) = range(A);
    
    % Duration
    duration(i) = stop - start + 1;
end

PVTimesDuration = peakVel .* duration ./ 1000;

%% PV vs. amplitude

fig2 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
tiledlayout(1, 1);
h = sgtitle(filename + ": Peak Velocity vs. Amplitude");
h.Interpreter = 'none';

ax4 = nexttile;
hold on;
scatter(amplitude, peakVel, 10, colormap(2, :), 'filled');
coefficients = polyfit(amplitude, peakVel, 1);
xFit = linspace(0, max(amplitude), 2);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'Color', colormap(2, :), 'HandleVisibility', 'off');
[r, pval] = corr(amplitude, peakVel);

grid on;
xlabel('Amplitude (deg)');
ylabel('PV (deg / s)');
xline(0, 'HandleVisibility', 'off');
yline(0, 'HandleVisibility', 'off');
legend(sprintf('y = %f * x + %f, R^2 = %f', coefficients(1), coefficients(2), r * r));

%% Duration vs. amplitude

fig3 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
tiledlayout(1, 1);
h = sgtitle(filename + ": Duration vs. Amplitude");
h.Interpreter = 'none';

ax3 = nexttile;
hold on;
scatter(amplitude, duration, 10, colormap(2, :), 'filled');
coefficients = polyfit(amplitude, duration, 1);
xFit = linspace(0, max(amplitude), 2);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'Color', colormap(2, :), 'HandleVisibility', 'off');
[r, pval] = corr(amplitude, duration);

grid on;
legend(sprintf('y = %f * x + %f, R^2 = %f', coefficients(1), coefficients(2), r * r));
xlabel('Amplitude (deg)');
ylabel('Duration (ms)');
xline(0, 'HandleVisibility', 'off');
yline(0, 'HandleVisibility', 'off');

%% PV * Duration vs. amplitude

fig4 = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
tiledlayout(1, 1);
h = sgtitle(filename + ": PV * Duration vs. Amplitude");
h.Interpreter = 'none';

ax5 = nexttile;
hold on;
scatter(amplitude, PVTimesDuration, 10, colormap(2, :), 'filled');
coefficients = polyfit(amplitude, PVTimesDuration, 1);
xFit = linspace(0, max(amplitude), 2);
yFit = polyval(coefficients, xFit);
plot(xFit, yFit, 'Color', colormap(2, :), 'HandleVisibility', 'off');
[r, pval] = corr(amplitude, PVTimesDuration);

grid on;
legend(sprintf('y = %f * x + %f, R^2 = %f', coefficients(1), coefficients(2), r * r));
xlabel('Amplitude (deg)');
ylabel('PV * Duration');
xline(0, 'HandleVisibility', 'off');
yline(0, 'HandleVisibility', 'off');

%% Save
saveas(fig2, strcat(filename, '_PV_vs_amplitude.fig'));
saveas(fig3, strcat(filename, '_duration_vs_amplitude.fig'));
saveas(fig4, strcat(filename, '_PV_times_duration_vs_amplitude.fig'));
save(strcat(filename, '_saccades.mat'), 'saccades', 'artifacts', 'isUsingCali', 'customScaleCh1', 'customScaleCh2', 'isRemovingSegment', 'removeTime');

if ~isUsingSaccades
    saveas(fig1, strcat(filename, '_traces.fig'));
end

end

%% Functions

% Remove transients
function datout = removeTransients(datin, thres)
    datout = datin;
    for i = 2 : length(datin) - 1
        if abs(datin(i) - datin(i - 1)) > thres && abs(datin(i) - datin(i + 1)) > thres && abs(datin(i - 1) - datin(i + 1)) < thres
            datout(i) = (datin(i - 1) + datin(i + 1)) / 2;
        end
    end
end

% Detect saccades
function saccadeCenters = detectSaccades(time, pos, vel, window, scaleFactor, isGraphing)
    [upper, lower] = getBounds(vel, window, scaleFactor);
    saccadeCenters = pos; % Copy of pos
    saccadeCenters(vel > lower & vel < upper) = NaN; % It's not a saccade if vel is within the bounds
    if isGraphing
        plot(time, upper, 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5); % Blue
        plot(time, lower, 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5);
    end
end

% Get start and stop
function datout = getStartStop(datin)
% The input should look like a vector (timeline) that has the saccade centers marked with the eye position
% The output is a matrix, where each row represents a saccade
%   Col 1 contains the start index of the saccade (inclusive)
%   Col 2 contains the end index of the saccade (inclusive)
    D = diff([false; ~isnan(datin); false]); % Marks start of saccade with 1, and end of saccade with -1
    starts = find(D > 0); % Col 1 contains the start index of each saccade
    stops = find(D < 0) - 1; % Col 2 contains the end index of each saccade
    datout = [starts, stops];
end

% Extend saccades to thres
function datout = extendSaccadesToThres(saccades, magnetPos, thres)
% The input should be a matrix, where each row represents a saccade
%   Col 1 contains the start index of the saccade
%   Col 2 contains the end index of the saccade
    datout = saccades;
    velocity = [diff(magnetPos); magnetPos(end)] * 1000; % Calculate the velocity without smoothing
    for i = 1:size(saccades, 1)
        start = saccades(i, 1);
        stop = saccades(i, 2);
        
        B = sign(velocity(round((start + stop) / 2))); % Sign of velocity at center of saccade
        
        % Extend to the left until velocity goes below thres
        A = velocity(1:start);
        if B > 0
            start = find(A <= thres, 1, 'last');
        else
            start = find(A >= -thres, 1, 'last');
        end
        if isempty(start)
            start = saccades(i, 1);
        end
        
        % Extend to the right until velocity goes below thres
        A = velocity(stop:end);
        if B > 0
            stop = find(A <= thres, 1, 'first') + stop - 1;
        else
            stop = find(A >= -thres, 1, 'first') + stop - 1;
        end
        if isempty(stop)
            stop = saccades(i, 2);
        end
        
        datout(i, 1) = start;
        datout(i, 2) = stop;
    end
end

% Put saccades in context
function datout = putInContext(saccades, magnetPos, merged)
    if merged
        % Returns a column
        % The column contains all the saccades, on the same time scale as magnetPos
        datout = NaN(size(magnetPos));
        for i = 1:size(saccades, 1)
            start = saccades(i, 1);
            stop = saccades(i, 2);
            datout(start:stop) = magnetPos(start:stop);
        end
    else
        % Returns a matrix
        % Each column contains one saccade
        % Each column is on the same time scale as magnetPos
        datout = NaN(numel(magnetPos), numel(saccades));
        for i = 1:size(saccades, 1)
            start = saccades(i, 1);
            stop = saccades(i, 2);
            datout(start:stop, i) = magnetPos(start:stop);
        end
    end
end

% Remove artifacts
function [saccades, artifacts] = removeArtifacts(saccades, magnetPos, maxDuration, neighborThres)
    remove = zeros(size(saccades, 1), 1);
    for i = 1:size(saccades, 1)
        % If it returns to starting pos within __ ms, then it's removed
        start = saccades(i, 1);
        stop = saccades(i, 2);
        if stop + maxDuration < numel(magnetPos)
            postSaccade = magnetPos(stop:stop + maxDuration);
        else
            postSaccade = magnetPos(stop:end);
        end
        if magnetPos(start) > min(postSaccade) && magnetPos(start) < max(postSaccade)
            remove(i) = 1;
        end

        % Consider the neighboring saccade before me
        % If I move back to near its starting position, then remove both of us
        % "Near" is within 1/2 of the amplitude of the neighbor
        if i >= 2
            displacement = abs(magnetPos(saccades(i - 1, 1)) - magnetPos(saccades(i, 2)));
            timeDiff = saccades(i, 1) - saccades(i - 1, 2); % Time from prev neighbor's end to my start
            
            prevAmp = range(magnetPos(saccades(i - 1, 1):saccades(i - 1, 2)));
            currAmp = range(magnetPos(saccades(i, 1):saccades(i, 2)));
            
            prevSign = sign(magnetPos(saccades(i - 1, 2)) - magnetPos(saccades(i - 1, 1)));
            currSign = sign(magnetPos(saccades(i, 2)) - magnetPos(saccades(i, 1)));
            
            if prevSign ~= currSign && timeDiff < neighborThres
                if displacement < 0.5 * prevAmp
                    remove(i) = 1;
                    remove(i - 1) = 1;
                elseif(currAmp < prevAmp)
                    remove(i) = 1;
                else
                    remove(i - 1) = 1;
                end
            end
        end
    end
    remove = logical(remove);
    
    artifacts = saccades(remove, :); % Returns the things that are removed
    saccades(remove, :) = [];
end

% Get bounds
function [upper, lower] = getBounds(vel, window, scaleFactor)
% Pass window=-1 in order to have a constant threshold
    if window > 0
        mads = movmad(vel, window);
        medians = movmedian(vel, window);
        upper = medians + scaleFactor * mads;
        lower = medians - scaleFactor * mads;
    else
        velMedian = median(vel);
        velMAD = mad(vel);
        upper = velMedian + velMAD * scaleFactor + zeros(size(vel));
        lower = velMedian - velMAD * scaleFactor + zeros(size(vel));
    end
end