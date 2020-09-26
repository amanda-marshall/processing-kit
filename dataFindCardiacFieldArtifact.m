function dataFindCardiacFieldArtifact(filePathECG, filePathComponents)

%specify path to fieldtrip on your drive
ft_defaults()

% Using ECG and Component data where noisy trials are already removed
% (so heart beat and CFA is easier to identify)
% Note: Actual removal of components needs to be done on original data

% filePathECG = '/Users/Amanda/Desktop/bioFeedback2/processed_data_ECG';
% filePathComponents = '/Users/Amanda/Desktop/bioFeedback2/processed_data_components';
outputPath = fileparts(filePathComponents);

% Load data
load(filePathECG,'data_ECG');
load(filePathComponents,'data_components');

% Identify heart beats relative to event of interest
% r-syncing makes identification of CFA easier - actual data is saved with
% original timing
sampleFrequency = data_ECG.fsample;
trialCount = length(data_components.trial);    % Number of trials
channelCount = length(data_components.label);  % Number of channels

% TODO: Why these specific numbers?
prePeakTime = 0.2;	% Time before each peak in seconds
postPeakTime = 0.6;	% Time after each peak in seconds

% Convert time (in sec) to samples, based on input sample frequency
prePeakSamples = prePeakTime * sampleFrequency;     % Samples before peak
postPeakSamples = postPeakTime * sampleFrequency;   % Samples after peak
totalSamples = prePeakSamples + 1 + postPeakSamples;      % Total sampling period

% Create empty array that can be used to mark trials that should be
% disregarded (0 == keep, 1 == exclude)
% exclude = zeros(trialCount, 1);

% TODO: Why these specific numbers?
% Note: Within each trial, this seems to be the region of interest. I.e.
% Jakob only considers hearbeats that happen between 1.5sec and 3.2sec in
% the trial. I am not sure why?
trialWindowOfInterest = [1.5 3.2]; % Time region of interest used for CFA identification

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
    rPeakTimes = ecg_timing(rWaveIndex);
    
    % Which of the R-Peaks fall within the window of interest?
    peakMatches = rPeakTimes > trialWindowOfInterest(1) & rPeakTimes < trialWindowOfInterest(2); % Creates an array of zeros (0) and ones (1) indicating which peaks are within the window of interest
    peakTimes = rPeakTimes(peakMatches); % Gets the timestamps representing the peaks of interest found above
    
    % If no peaks are identified within the region of interest?
    if isempty(peakTimes)
        
        % Print warning
        fprintf('Warning: No hearbeats found in trial %d\n', trialNumber);
        
        %{
        % Mark trial for exclusion from analysis
        exclude(trialNumber) = 1;
        
        % TODO: This looks like Jakob is writing zeros to the "bad" trials
        % in the ECG and Components dataset. However, he seems to
        % completely remove them after the loop anyway. Is this needed?
        data_ECG.trial{trialNumber} = zeros(1, totalSamples);
        data_components.trial{trialNumber} = zeros(channelCount, totalSamples);
        %}
       
    % Else, if there are peaks within the region of interest?
    else
        
        % Get time of first peak within region of interest
        firstPeakTime = peakTimes(1);
        
        % If the end of the first peak period is longer than the trial, the
        % trial is invalid?
        if firstPeakTime + postPeakTime > ecg_timing(end)
            
            % Print warning
            fprintf('Warning: Heartbeat in trial %d is too late for extraction\n', trialNumber);
            
            %{
            % Mark this trial for exclusion from analysis
            exclude(trialNumber) = 1;
            
            % TODO: This looks like Jakob is writing zeros to the "bad" trials
            % in the ECG and Components dataset. However, he seems to
            % completely remove them after the loop anyway. Is this needed?
            data_ECG.trial{trialNumber} = zeros(1, totalSamples);
            data_components.trial{trialNumber} = zeros(channelCount, totalSamples);
            %}
            
        else
            
            % Get the sample indeces representing R-Peaks
            peakSampleIndeces = rWaveIndex(peakMatches);
            
            for sampleIndex = 1 : length(peakSampleIndeces)
                
                % Define the start sample and end sample around the first
                % R-Peak within the window of interest, based on the
                % prePeakTime and postPeakTime defined at the top
                peakStartSampleIndex = peakSampleIndeces(sampleIndex) - (prePeakSamples + 1);
                peakEndSampleIndex = peakSampleIndeces(sampleIndex) + postPeakSamples - 1;

                % Replace trial data with ECG data around the first peak within
                % the window of interest
                trialECGData = data_ECG.trial{trialNumber}(:, peakStartSampleIndex:peakEndSampleIndex);
                trialComponentData = data_components.trial{trialNumber}(:, peakStartSampleIndex:peakEndSampleIndex);

                trials_ECG{end + 1} = trialECGData;
                trials_components{end + 1} = trialComponentData;
                
            end;
            
        end
    
    end
    
    %{
    % Generate new equidistant time stamps for samples
    sampleTimes = linspace(-1 * prePeakTime, postPeakTime, totalSamples);
        
    % Generate timings for the ECG trial, based on the prePeakTime and
    % postPeakTime defined at the top
    data_ECG.time{trialNumber} =  sampleTimes;
    data_ECG.sampleinfo(trialNumber, :) = [1 totalSamples]; % TODO: Why this specific sample info?
    
    % Generate timings for the Components trial, based on the prePeakTime
    % and postPeakTime defined at the top
    data_components.time{trialNumber} = sampleTimes;
    data_components.sampleinfo(trialNumber, :) = [1 totalSamples]; % Todo: Why this specific sample info?
    %}
    
    
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
for cnum=4:-1:1
    
    componentId = componentIdsSorted(cnum);
    otherComponentIds = setdiff(componentIdsSorted, componentId);
    
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
    
    filePath = fullfile(outputPath, ['component_', num2str(componentId)]);
    saveas(fig, filePath, 'png');
     
end

disp('****************************************************************************************')


