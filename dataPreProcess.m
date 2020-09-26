clear all
plist = [12];
addpath('F:\oneDrive\lmuprojects\fieldtrip-20181029')
ft_defaults()
for i=1:length(plist)
    clear eegdata components exclIndex
    p = plist(i);
    pnum = num2str(p);
    if p < 10
        pnum = ['0' pnum];
    end
    
    disp(['participant ' pnum ': importing to fieldtrip--------------------']);
    cfg = [];
    cfg.datafile     = ['eeg1_raweeg/motorLearn00' pnum '.eeg'];
    dat_raw = ft_preprocessing(cfg);
%     
%     %this label was spelled incorrectly in first data sets
%     dat_raw.label{1} = 'FP1';


  
    if exist('exclIndex','var')
        numchannels =length(dat_raw.label);
        cfg         = [];
        cfg.channel = setdiff(1:numchannels, exclIndex);
        dat_raw        = ft_selectdata(cfg, dat_raw);
    end

%adding the implicit reference channel to the data
   
    %last channel = ECG
    numchannels =length(dat_raw.label);
    cfg         = [];
    cfg.channel = 1:numchannels-1;
    dataNoECG        = ft_selectdata(cfg, dat_raw);
    %switching to average reference of all channels
    cfg = [];
    cfg.implicitref   = 'FCz';
    cfg.reref         = 'yes';
    cfg.refchannel    = 'all';
    cfg.refmethod     = 'avg';
    dataNoECG = ft_preprocessing(cfg,dataNoECG);
    
    
    

    
    
    
    cfg         = [];
    cfg.channel = {'ECG'};
    ecgOnly      = ft_selectdata(cfg, dat_raw);
    dat_raw =      ft_appenddata(cfg, dataNoECG,ecgOnly);
    
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 40;
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    dat_raw = ft_preprocessing(cfg,dat_raw);
    
    
    
    
    cfg = [];
    cfg.dataset             = ['eeg1_raweeg/motorLearn00' pnum '.vhdr'];
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = {'S111','S112','S113','S114','S121','S122','S123','S124'};
    cfg.trialdef.prestim    = 1.5;
    cfg.trialdef.poststim    = 4.5;
    trialinfo        = ft_definetrial(cfg);
    trl = trialinfo.trl;
    
    eegdata = ft_redefinetrial(trialinfo, dat_raw);
    
    
    %downsampling and switching to single precision since this is more
    %economically
    cfg = [];
    cfg.resamplefs = 250;
    cfg.detrend    = 'no';
    eegdata = ft_resampledata(cfg, eegdata);
    %     smplpos = 1;
    numtrials = size(eegdata.trial,2);
    for tnum=1:numtrials
        %         trl(tnum,1) = smplpos;
        trl(tnum,2) = trl(tnum,1) + size(eegdata.trial{tnum},2)-1;
        %         smplpos = trl(tnum,2)+1;
        trl(tnum,3) = -1*(find(eegdata.time{1}>=0,1)-1); %all samples before point 0
    end
    % eegdata.cfg.trl = trl;
    cfg = [];
    cfg.precision = 'single';
    eegdata = ft_preprocessing(cfg,eegdata);
    
    %
    clear dataNoECG ecgOnly dat_raw
    
    %no splitting up data again, because ECG is not part of the ICA
        %last channel = ECG
    numchannels =length(eegdata.label);
    cfg         = [];
    cfg.channel = 1:numchannels-1;
    dataNoECG        = ft_selectdata(cfg, eegdata);
        
    cfg         = [];
    cfg.channel = {'ECG'};
    ecgOnly      = ft_selectdata(cfg, eegdata);
    
    
    cfg        = [];
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
    cfg.runica.extended = 0;
    components = ft_componentanalysis(cfg, dataNoECG);
    
    eegdata = dataNoECG;
    eegdata.cfg.trl = trl;
    save(['eeg2_ICA/learnEEG' num2str(p)], 'eegdata');
   
    save(['eeg2_ICA/learnCOMP' num2str(p)], 'components');
     save(['eeg2_ICA/learnECG' num2str(p)], 'ecgOnly');
    
end