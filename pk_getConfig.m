function [experimentConfig, participantConfig] = pk_getConfig(configFilePath, participantId)

  if ~exist(configFilePath, 'file')
    error('Cannot open configuration. File %s could not be found.', configFilePath);
  end;

  experimentConfig = jsondecode(fileread(configFilePath));
  participantConfig = [];

  % If config could be loaded and a participant ID was specified…
  if ~isempty(experimentConfig) && ~isempty(participantId)
    pKey = strcat('p', num2str(participantId));
    if isfield(experimentConfig, 'participants') && isfield(experimentConfig.participants, pKey)
      participantConfig = experimentConfig.participants.(pKey);
    end;
  end;

end
