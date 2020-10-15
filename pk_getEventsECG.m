function events = pk_getEventsECG(dataECG)

  sampleRate = dataECG.fsample;
  trials = dataECG.trial;
  trialCount = length(trials);

  events = struct('type', {}, 'value', {}, 'sample', {}, 'duration', {}, 'offset', {});

  % Loop over trials in dataset
  for trialIndex = 1 : trialCount

    % Get ECG data + sample times for current trial
    trialSampleBounds = dataECG.sampleinfo(trialIndex, :);
    trialSampleStart = trialSampleBounds(1);
    trialTimes = dataECG.time{trialIndex};
    trialData = dataECG.trial{trialIndex}; % ECG data for this trial
    trialData = -1 * trialData;

    % Detect R-peaks in ECG trace for current trial
    [~, peakSamplesIndeces] = pan_tompkin(trialData, sampleRate, 0);

    for i = 1 : length(peakSamplesIndeces)

      % Get R-peak sample
      sampleIndex = peakSamplesIndeces(i);
      sampleTime = trialTimes(sampleIndex);

      sample = trialSampleStart + sampleTime * sampleRate;

      % Create event
      eventTemplate = [];
      eventTemplate.type = 'Stimulus';
      eventTemplate.value = 'R-Peak';
      eventTemplate.sample = sample;
      eventTemplate.duration = 1;
      eventTemplate.offset = [];

      % Add event to list
      events(end+1) = eventTemplate;

    end

  end

  fprintf('Found %d R-Peak events in the ECG trace\n', length(events));

end
