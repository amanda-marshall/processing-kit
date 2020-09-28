function pk_dataPreProcess(filePath, cfgSegmentation)

outputPath = fileparts(filePath);

% Settings
prefs = [];
prefs.eegFile = filePath; % File to read
prefs.ecgChannel = 'ECG';
prefs.eegChannels = { 'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6', 'T7', 'C3', 'Cz', 'C4', 'T8', 'TP9', 'CP5', 'CP1', 'CP2', 'CP6', 'TP10', 'P7', 'P3', 'Pz', 'P4', 'P8', 'PO9', 'O1', 'Oz', 'O2', 'PO10', 'AF7', 'AF3', 'AF4', 'AF8', 'F5', 'F1', 'F2', 'F6', 'F9', 'FT7', 'FC3', 'FC4', 'FT8', 'F10', 'C5', 'C1', 'C2', 'C6', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'P5', 'P1', 'P2', 'P6', 'PO7', 'POz', 'PO3', 'P03', 'PO4', 'PO8' };
% NOTE: This channel configuration contains both the correct PO3 and
% incorrect P03 channel (P ZERO 3), for misconfigured datasets.

% Read raw EEG data into MATLAB
cfg = [];
cfg.dataset = prefs.eegFile;
cfg.channel = prefs.eegChannels;

% Read in EEG Data
fprintf('\n### 1. Reading EEG Data\n\n');
data_EEG = ft_preprocessing(cfg);

% Read in ECG Data
fprintf('\n### 2. Reading ECG Data\n\n');
cfg.channel = prefs.ecgChannel;
data_ECG = ft_preprocessing(cfg);

% Set up EEG filters
cfg = [];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfreq = 1;
cfg.lpfreq = 40;

% Set up EEG re-referencing
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';

% Filter and re-reference EEG Data
fprintf('\n### 3. Filtering + Re-referencing EEG Data\n\n');
data_EEG = ft_preprocessing(cfg, data_EEG);

% NOTE: We used to segment the data at this point, but we will instead do this
% after filtering, re-referencing, component analysing and CFA detection.

% Combine EEG and ECG Datasets
fprintf('\n### 4. Combine EEG + ECG Datasets\n\n');
cfg = [];
data = ft_appenddata(cfg, data_EEG, data_ECG);

% Create a trial definition configuration
fprintf('\n### 5. Define trial windows for segmentation\n\n');
cfg = [];
cfg.dataset = prefs.eegFile;
cfg.trialfun = 'ft_trialfun_bracketed'; % Custom trial definition function that uses triggers as markers for the start and end of the trial
cfg.trialdef.eventtype = 'Stimulus'; % Trial event type
cfg = pk_mergeStructs(cfg, cfgSegmentation); % Merge passed trial definition with config object
trial_definition = ft_definetrial(cfg);

% Segment data based on trial definition
fprintf('\n### 6. Segment merged dataset into full trial windows\n\n');
data_segmented = ft_redefinetrial(trial_definition, data);

% Split data into EEG and ECG after segmentation
fprintf('\n### 7. Split dataset into EEG and ECG after segmentation\n\n');
cfg = [];
cfg.channel = prefs.eegChannels;
data_EEG = ft_selectdata(cfg, data_segmented);
cfg.channel = prefs.ecgChannel;
data_ECG = ft_selectdata(cfg, data_segmented);

% Run Independent Component Analysis
fprintf('\n### 8. Run Independent Component Analysis\n\n');
cfg = [];
cfg.method = 'runica'; % This is the default and uses the implementation from EEGLAB
cfg.runica.extended = 0;
data_components = ft_componentanalysis(cfg, data_EEG);

% Save processed data
fprintf('\n### 9. Save data to disk\n\n');
filePathEEG = fullfile(outputPath, 'processed_data_EEG');
filePathECG = fullfile(outputPath, 'processed_data_ECG');
filePathComponents = fullfile(outputPath, 'processed_data_components');
save(filePathEEG, 'data_EEG');
save(filePathECG, 'data_ECG');
save(filePathComponents, 'data_components');

% Done
fprintf('\n### Done\n\n');
