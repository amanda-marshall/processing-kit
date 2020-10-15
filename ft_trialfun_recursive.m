function [trl, events] = ft_trialfun_recursive(cfg)

  % this function requires the following fields to be specified
  % cfg.dataset
  % cfg.trialdef.eventtype
  % cfg.trialdef.eventvalue
  % cfg.trialdef.prestim
  % cfg.trialdef.poststim

  hdr   = ft_read_header(cfg.dataset);
  events = ft_read_event(cfg.dataset);

  % If a custom events struct array was supplied, merge it with the ones
  % from the EEG file and sort them according to their samples
  hasCustomEvents = isfield(cfg, 'events');
  if hasCustomEvents
    events = [events, cfg.events];
    eventsTable = struct2table(events);
    eventsTableSorted = sortrows(eventsTable, 'sample');
    events = table2struct(eventsTableSorted);
  end

  cfg.sampleRate = hdr.Fs;
  cfg.events = events;
  cfg.indexStart = 1;
  cfg.indexEnd = length(events);

  trl = [];
  trl = pk_find_trials(cfg, trl);

end
