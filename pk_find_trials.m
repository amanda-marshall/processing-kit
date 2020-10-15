function trl = pk_find_trials(cfg, passedTrl)

  % Deconstruct passed configuration
  trl = passedTrl;
  trialDefinition = cfg.trialdef;
  sampleRate = cfg.sampleRate;
  events = cfg.events;
  indexStart = cfg.indexStart;
  indexEnd = cfg.indexEnd;

  hasEndEvent = isfield(trialDefinition, 'eventvalueend');
  isNestedTrial = isfield(trialDefinition, 'trialdef');
  hasNthTarget = isfield(trialDefinition, 'nthTarget');

  nthCount = 0;

  % Loop over all events in the current subset
  for indexStartEvent = indexStart:indexEnd
    startEvent = events(indexStartEvent);

    % Test if event is a valid start event
    isStartEventTypeMatch = strcmp(startEvent.type, trialDefinition.eventtype);
    isStartEventValueMatch = isStartEventTypeMatch && ismember(startEvent.value, trialDefinition.eventvalue);
    isStartEvent = isStartEventTypeMatch && isStartEventValueMatch;

    % Continue if event matches type and value of a start event
    if isStartEvent

      nthCount = nthCount + 1;

      % Abort loop if the nth target has been found
      if hasNthTarget && nthCount > trialDefinition.nthTarget
        break;
      end

      % Variables to record for each trial
      sampleStart = startEvent.sample - trialDefinition.prestim * sampleRate;
      sampleEnd = max(sampleStart, sampleStart + trialDefinition.poststim * sampleRate - 1);
      offset = sampleStart - startEvent.sample;

      % fprintf('Sample start is %d, sampleEnd is %d, offset is %d\n', sampleStart, sampleEnd, offset);

      if hasEndEvent && indexStartEvent < indexEnd
        for indexEndEvent = indexStartEvent + 1 : indexEnd

          endEvent = events(indexEndEvent);

          % Test if event is a valid end event
          isEndEventTypeMatch = strcmp(endEvent.type, trialDefinition.eventtype);
          isEndEventValueMatch = isEndEventTypeMatch && ismember(endEvent.value, trialDefinition.eventvalueend);
          isEndEvent = isEndEventTypeMatch && isEndEventValueMatch;

          if isEndEvent
            sampleEnd = endEvent.sample + trialDefinition.poststim * sampleRate - 1;
            break;
          end

        end
      end

      % If the trial definition isn't nested, record the trial info and move on
      if ~isNestedTrial

        % Skip this trial, if it doesn't meet the nth target requirement
        if hasNthTarget && trialDefinition.nthTarget ~= nthCount
          continue;
        end

        % Add trial data to output
        trl(end+1, :) = [round([sampleStart sampleEnd offset])];
        continue;

      end;

      % Find earliest event index in trial definition
      nestedStartIndex = indexStartEvent;
      for indexReverse = 1:indexStartEvent
        index = indexStartEvent - (indexReverse - 1);
        testEvent = events(index);
        if testEvent.sample >= sampleStart
          nestedStartIndex = index;
        else
          break;
        end
      end

      % Find latest event index in trial definition
      nestedEndIndex = indexStartEvent;
      for index = indexStartEvent:indexEnd
        testEvent = events(index);
        if testEvent.sample <= sampleEnd
          nestedEndIndex = index;
        else
          break;
        end
      end

      % Recursively resolve nested trials
      cfg.trialdef = trialDefinition.trialdef;
      cfg.indexStart = nestedStartIndex;
      cfg.indexEnd = nestedEndIndex;
      trl = pk_find_trials(cfg, trl);

    end

  end

end
