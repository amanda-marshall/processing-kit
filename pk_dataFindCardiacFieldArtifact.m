function pk_dataFindCardiacFieldArtifact(filePathECG, filePathComponents, passedPrefs)

% Internal preferences
prefs = [];
prefs.shouldSkipPreviouslyProcessed = false;
prefs.outputComponentCount = 4; % Home many components should be written out?
prefs.outputPath = fileparts(filePathComponents);
prefs.prePeakTime = 0.2;
prefs.postPeakTime = 0.6;
prefs = pk_mergeStructs(prefs, passedPrefs); % Merge passed prefs with defaults

existCount = 0;
outputFileNames = cell(1, prefs.outputComponentCount);
for i = 1:prefs.outputComponentCount
  figureFilePath = fullfile(prefs.outputPath, ['component_', num2str(i), '.png']);
  if exist(figureFilePath, 'file')
    existCount = existCount + 1;
  end
  outputFileNames(i) = {figureFilePath};
end

% Bail, if all output files already exist
if prefs.shouldSkipPreviouslyProcessed && existCount == prefs.outputComponentCount
  fprintf('\n### Output Data for %s already exists, skipping\n\n', filePathComponents);
  return;
end;

% Using ECG and Component data where noisy trials are already removed
% (so heart beat and CFA is easier to identify)
% Note: Actual removal of components needs to be done on original data

% Load data
load(filePathECG, 'data_ECG');
load(filePathComponents, 'data_components');

% Identify heart beats relative to event of interest
% r-syncing makes identification of CFA easier - actual data is saved with
% original timing
sampleFrequency = data_ECG.fsample;
trialCount = length(data_components.trial);    % Number of trials

prePeakTime = prefs.prePeakTime;	% Time before each peak in seconds
postPeakTime = prefs.postPeakTime;	% Time after each peak in seconds

% Convert time (in sec) to samples, based on input sample frequency
prePeakSamples = prePeakTime * sampleFrequency;     % Samples before peak
postPeakSamples = postPeakTime * sampleFrequency;   % Samples after peak
totalSamples = prePeakSamples + 1 + postPeakSamples;      % Total sampling period

trials_ECG = {};
trials_components = {};

% Loop over individual trials
for trialNumber = 1 : trialCount

    ecg_timing = data_ECG.time{trialNumber}; % Array of time stamps for each sample in the trial
    ecg_trial = data_ECG.trial{trialNumber}; % ECG data for this trial
    ecg_trial = -1 * ecg_trial; % TODO: Why is this multiplied by -1?

    % Run the Pan Tompkins algorithm over the ECG data for this trial to
    % detect R peaks (listed in the returned index)
    [~, rWaveIndex] = pan_tompkin(ecg_trial, sampleFrequency, 0);

    % At which time (in s) do R-Peaks appear in this trial?
    peakTimes = ecg_timing(rWaveIndex);

    % If no peaks are identified within the region of interest?
    if isempty(peakTimes)

        % Print warning
        fprintf('Warning: No hearbeats found in trial %d\n', trialNumber);

    % Else, if there are peaks within the region of interest?
    else

        % Get the sample indeces representing R-Peaks
        peakSampleIndeces = rWaveIndex;

        for sampleIndex = 1 : length(peakSampleIndeces)

            % Define the start sample and end sample around the first
            % R-Peak within the window of interest, based on the
            % prePeakTime and postPeakTime defined at the top
            peakStartSampleIndex = peakSampleIndeces(sampleIndex) - (prePeakSamples + 1);
            peakEndSampleIndex = peakSampleIndeces(sampleIndex) + postPeakSamples - 1;

            % Skip this heartbeat, if the sample window too early/late
            if peakStartSampleIndex < 1 || peakEndSampleIndex > length(data_ECG.trial{trialNumber})
                continue;
            end;

            % Replace trial data with ECG data around the first peak within
            % the window of interest
            trialECGData = data_ECG.trial{trialNumber}(:, peakStartSampleIndex:peakEndSampleIndex);
            trialComponentData = data_components.trial{trialNumber}(:, peakStartSampleIndex:peakEndSampleIndex);

            trials_ECG{end + 1} = trialECGData;
            trials_components{end + 1} = trialComponentData;

        end;

    end

end

% Replace original trials with new ones
data_ECG.trial = trials_ECG;
data_components.trial = trials_components;

% How many trials did we end up generating?
trialCount = length(trials_ECG);

% Generate new equidistant time stamps for samples
sampleTrialTimes = linspace(-1 * prePeakTime, postPeakTime, totalSamples);
sampleTrialTimesMatrix = repmat({sampleTrialTimes}, 1, trialCount);
data_ECG.time = sampleTrialTimesMatrix;
data_components.time = sampleTrialTimesMatrix;

% Generate new sample info
sampleInfoMatrix = repmat([1 totalSamples], trialCount, 1);
data_ECG.sampleinfo = sampleInfoMatrix;
data_components.sampleinfo = sampleInfoMatrix;

% Exclude all trials where no suitable peak could be found
% Afterward, combine ECG + components dataset into one dataset
cfg = [];
cfg.trials = linspace(1, trialCount, trialCount); % setdiff(1:trialCount, find(exclude));
data_rSync_ECG = ft_selectdata(cfg, data_ECG);
data_rSync_components = ft_selectdata(cfg, data_components);
data_rSync = ft_appenddata([], data_rSync_components, data_rSync_ECG);

% Compute a frequency decomposition of all components and the ECG
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.foilim     = [1 40];
cfg.taper      = 'hanning';
cfg.pad        = 'maxperlen';
freq           = ft_freqanalysis(cfg, data_rSync);

% Compute coherence between all components and the ECG
cfg            = [];
cfg.channelcmb = {'all' 'ECG'};
cfg.jackknife  = 'no';
cfg.method     = 'coh';
fdcomp         = ft_connectivityanalysis(cfg, freq);

[~, componentIdsSorted] = sort(mean(fdcomp.cohspctrm, 2), 'descend');

% How many components do we want to see?
componentIdsOfInterest = componentIdsSorted(1:prefs.outputComponentCount);
for cnum = prefs.outputComponentCount:-1:1

    componentId = componentIdsSorted(cnum);
    otherComponentIds = setdiff(componentIdsOfInterest, componentId);

    fig = figure('visible', 'off');

    subplot(2,3,1);
    cfg = [];
    cfg.component = componentId;       % specify the component(s) that should be plotted
    cfg.layout    = 'layout65.mat'; % specify the layout file that should be used for plotting
    ft_topoplotIC(cfg, data_components)
    subplot(2,3,4);
    hold on
    plot(fdcomp.freq, abs(fdcomp.cohspctrm(componentId,:)));
    plot(fdcomp.freq, abs(fdcomp.cohspctrm(otherComponentIds,:)),':');
    title(['comp-ECG coherence: ' num2str(round(mean(fdcomp.cohspctrm(componentId,:),2),2))])
    hold off
    xlabel('frequency')
    ylabel('coherence')

    cfg           = [];
    timelock      = ft_timelockanalysis(cfg, data_rSync);
    cfg = [];
    cfg.toilim    =  [-0.2 0.6];
    timelock = ft_redefinetrial(cfg, timelock);
    subplot(2,3,2);
    plot(timelock.time, timelock.avg(end,:)) %plotting ECG
    title('average ECG')
    xlabel('time')
    ylabel('activity')

    subplot(2,3,5);
    hold on
    plot(timelock.time, timelock.avg(componentId,:));
    plot(timelock.time, timelock.avg(otherComponentIds,:),':')
    hold off
    title('average component')
    xlabel('time')
    ylabel('activity')


    cfg = [];
    cfg.toilim    =  [-0.2 0.6];
    componentsCut = ft_redefinetrial(cfg, data_rSync_components);
    ecgCut = ft_redefinetrial(cfg, data_rSync_ECG);

    clear tdat cdat
    for trialNumber=1:length(componentsCut.trial)
        tdat(trialNumber,:) = componentsCut.trial{trialNumber}(componentId,:);
        cdat(trialNumber,:) = ecgCut.trial{trialNumber};
    end
    subplot(2,3,3); imagesc(componentsCut.time{1},1:length(componentsCut.trial),cdat);
    title('trial-wise ECG activity')
    xlabel('time')
    ylabel('trial')

    subplot(2,3,6); imagesc(componentsCut.time{1},1:length(componentsCut.trial),tdat);
    title('trial-wise comp activity')
    xlabel('time')
    ylabel('trial')

    saveas(fig, char(outputFileNames(cnum)));

end

disp('****************************************************************************************')
