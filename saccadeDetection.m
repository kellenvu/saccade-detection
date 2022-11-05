%{
Code by Kellen Vu.

(1) Put a .smr file in this directory that contains eye position data.

(2) Put a cali.mat file in this directory that contains the scale factor
for the eye position data.

(3) Run this program. It will use a thresholding algorithm to detect the
saccades, and then spit out a data.mat file that contains an X vector (eye
position) and a Y vector (labels for saccade vs. non-saccade).
%}

tic

%% Options

isUsingSaccades = false;

isUsingCali = true;
customScaleCh1 = 95.05;
customScaleCh2 = 0;

% Manually exclude time segments from saccade detection
isExcludingTime = false;
excludeTime = [ % Each row is (start, stop) in seconds
    700, 1470;
];

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

% Check for .smr file
if isempty(dir('*.smr'))
    error('Make sure this .m file is in the same directory as the .smr file')
end

% Load mat files
if isUsingSaccades
    load(saccadeFile.name)
end
if isUsingCali && ~exist('caliFile', 'var')
    error('Make sure this .m file is in the same directory as the cali.mat file')
end
if isUsingCali
    load(caliFile.name, 'scaleCh1', 'scaleCh2')
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
tiledlayout(2, 1)
h = sgtitle(filename);
h.Interpreter = 'none';
colormap = brewermap(4, 'Set1');

% Timescale
magnetTime = (0:0.001:length(magnetPos) / 1000 - 0.001)';

% Magnet pos
ax1_1 = nexttile;
hold on
plot(magnetTime, magnetPos, 'Color', [0.7, 0.7, 0.7]);

ylabel('Position (deg)')
title('Position Trace')

% Magnet vel
ax1_2 = nexttile;
hold on
magnetVel = movingslope(magnetPos, 30) * 1000;
plot(magnetTime, magnetVel, 'Color', [0.7, 0.7, 0.7]);

% Detect saccades
if ~isUsingSaccades
    % Saccades
    saccadesInContext = detectSaccades(magnetTime, magnetPos, magnetVel, 1000, 5, false);
    saccades = getStartStop(saccadesInContext);
    
    if isExcludingTime
        for i = 1:size(excludeTime, 1)
            start = excludeTime(i, 1);
            stop = excludeTime(i, 2);
            removeIndices = saccades(:, 1) > start * 1000 & saccades(:, 2) < stop * 1000;
            saccades(removeIndices, :) = [];
        end
    end
    
    % Get peak vel
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
    [saccades, artifacts] = removeArtifacts(saccades, magnetPos, 100, 50);
end

ylabel('Velocity (deg / s)')
title('Velocity')
legend('Velocity')

% Plot saccades
axes(ax1_1)
saccadesInContext = putInContext(saccades, magnetPos, true);
plot(magnetTime, saccadesInContext, 'Color', colormap(2, :), 'LineWidth', 1)
artifactsInContext = putInContext(artifacts, magnetPos, true);
plot(magnetTime, artifactsInContext, 'Color', colormap(1, :), 'LineWidth', 1)

grid on
xlabel('Time (s)')
ylabel('Position (deg)')
title('Position')
legend('Pos', 'Saccades', 'Non-saccades')

linkaxes([ax1_1, ax1_2], 'x')

%% Save

save(filename + "_saccades.mat", 'saccades', 'artifacts', 'isUsingCali', 'customScaleCh1', 'customScaleCh2', 'isExcludingTime', 'excludeTime')

beep
toc

%% Functions

function datout = removeTransients(datin, thres)
    % Detect every transient (an outlier that is exactly one data point), and replace it with the average of its two neighbors.
    % :param datin: A vector of data over time
    % :param thres: The absolute minimum threshold for the amplitude of a transient
    % :return datout: The data with the transients removed
    datout = datin;
    for i = 2:length(datin) - 1
        if abs(datin(i) - datin(i - 1)) > thres && abs(datin(i) - datin(i + 1)) > thres && abs(datin(i - 1) - datin(i + 1)) < thres
            datout(i) = (datin(i - 1) + datin(i + 1)) / 2;
        end
    end
end

function saccadesInContext = detectSaccades(time, pos, vel, window, scaleFactor, isGraphing)
    % Detect saccades, and return them in a plottable vector.
    % :param time: A vector of time data
    % :param pos: A vector of position data
    % :param vel: A vector of velocity data
    % :param window: The size of the sliding window (ms) for determining the velocity threshold of a sacccade
    % :param scaleFactor: The scale factor for determining the velocity threshold of a saccade
    % :param isGraphing: If true, this plots the velocity threshold on the current axes
    % :return saccadesInContext: A plottable vector of the saccades, where all the non-saccade data are NaN
    [upper, lower] = getBounds(vel, window, scaleFactor);
    saccadesInContext = pos; % Copy of pos
    saccadesInContext(vel > lower & vel < upper) = NaN; % It's not a saccade if vel is within the bounds
    if isGraphing
        plot(time, upper, 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5) % Blue
        plot(time, lower, 'Color', [0 0.4470 0.7410], 'LineWidth', 0.5)
    end
end

function datout = getStartStop(datin)
    % Given some events "in context," find the start and stop index of each event.
    % :param datin: A vector of event data in context (i.e. all the non-event data are NaN)
    % :return datout: A matrix, where each row represents a saccade
    %   Col 1 contains the start index of the saccade (inclusive)
    %   Col 2 contains the end index of the saccade (inclusive)
    D = diff([false; ~isnan(datin); false]); % Marks start of saccade with 1, and end of saccade with -1
    starts = find(D > 0); % Col 1 contains the start index of each saccade
    stops = find(D < 0) - 1; % Col 2 contains the end index of each saccade
    datout = [starts, stops];
end

function datout = extendSaccadesToThres(saccades, pos, thres)
    % Extend the start and stop indices of every saccade.
    % :param saccades: Matrix of start and stop indices for saccades
    % :param pos: A vector of position data
    % :param thres: The absolute maximum velocity threshold for what should be considered a saccade
    % :return datout: Matrix of start and stop indices, extended
    datout = saccades;
    velocity = [diff(pos); pos(end)] * 1000; % Calculate the velocity without smoothing
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
    datout = unique(datout, 'rows');
end

function datout = putInContext(saccades, pos, isMerged)
    % Given the start and stop indices for events, return a plottable vector of those events.
    % :param saccades: Matrix of start and stop indices for saccades
    % :param pos: Vector of position data
    % :isMerged: Boolean for whether the output should be merged into one vector
    % :return datout: A plottable vector of the events, where all the non-event data are NaN
    if isMerged
        % Return a vector
        datout = NaN(size(pos));
        for i = 1:size(saccades, 1)
            start = saccades(i, 1);
            stop = saccades(i, 2);
            datout(start:stop) = pos(start:stop);
        end
    else
        % Return a matrix
        % Each column contains one saccade
        datout = NaN(numel(pos), numel(saccades));
        for i = 1:size(saccades, 1)
            start = saccades(i, 1);
            stop = saccades(i, 2);
            datout(start:stop, i) = pos(start:stop);
        end
    end
end

% Remove artifacts
function [saccades, artifacts] = removeArtifacts(saccades, magnetPos, maxDuration, neighborThres)
    % Given the start and stop indices for saccades, remove the rows that are artifacts, and return those artifacts separately.
    % :param saccades: Matrix of start and stop indices for saccades
    % :param pos: Vector of position data
    % :param maxDuration: If it returns to the starting pos within __ ms, then it's removed
    %   Pass a negative number to make the maxDuration a multiple of the duration of the saccade
    % :param neighborThres: Absolute max threshold for how far to look for neighboring saccades (ms)
    % :return saccades: Matrix of start and stop indices for saccades, with artifacts removed
    % :return artifacts: Matrix of start and stop indices for artifacts
    remove = zeros(size(saccades, 1), 1);
    if maxDuration < 0
        maxDuration = (stop - start) * abs(maxDuration);
    end
    for i = 1:size(saccades, 1)
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
    % Compute the min velocity threshold for what should be considered a saccade, using a sliding window.
    % :param vel: A vector of velocity data
    % :param window: The size of the moving window (ms)
    %   Pass window=-1 to use a constant threshold
    % :param scaleFactor: The scale factor used to compute the threshold
    %   Higher means the threshold is higher
    % :return upper: A vector of the upper velocity threshold
    % :return lower: A vector of the lower velocity threshold
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