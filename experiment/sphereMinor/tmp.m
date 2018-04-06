cd('/Users/marcw/Dropbox/work/fNIRS Hackathon/data/artinis/Measurement Examples')

cfg = [];
cfg.dataset = '1102 24 channel_2.oxy3';

[data] = ft_preprocessing(cfg);
cfg = [];
ft_databrowser(cfg, data);

%%

cfg = [];
cfg.dataset = '1102 24 channel_2.oxy3';
cfg.trialdef = [];
cfg.trialdef.eventtype = '?';


ft_definetrial(cfg);



%%


cfg.trialdef.eventtype  = 'event';
cfg.trialdef.eventvalue = 'A';
cfg.trialdef.prestim    = 10;
cfg.trialdef.poststim   = 35;
cfg = ft_definetrial(cfg);
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.dpf = 5.9;
conc_data = ft_transform_ODs(cfg, data);

