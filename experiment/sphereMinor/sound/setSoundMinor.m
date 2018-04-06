function setSoundMinor(snd,cfg,RZ6str)
% SETSOUND(SND,CFG,RP2STR)
%
% Set all parameters for sound presentation

Z			= snd.Z;
warning('Max sound level needs to be determined');
atten		= max(90-snd.intensity,0);
dur			= snd.offdelay-snd.ondelay;
sndsetup	= cfg.lookup(Z+1,2:3);

%% 
% cfg.(RZ6str).SetTagVal('delaySND1',snd.ondelay);  %
% cfg.(RZ6str).SetTagVal('delaySND2',snd.ondelay);  %
% cfg.(RZ6str).SetTagVal('soundDur1',dur);  % default?
% cfg.(RZ6str).SetTagVal('soundDur2',dur);  % default?

cfg.RZ6_1.SetTagVal('delaySND',snd.ondelay);
cfg.RZ6_1.SetTagVal('soundDur',dur);  %

cfg.RZ6_1.SetTagVal('AttenuationA',atten);
cfg.RZ6_1.SetTagVal('AttenuationB',atten);

cfg.RZ6_1.SetTagVal('eventSND',snd.onevent+1);

MUX(cfg.(RZ6str),sndsetup(1),sndsetup(2));
