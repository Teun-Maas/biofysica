function [RT, TN, D, A, RV, RD, MD, RW] = pa_spk_readdat(fname)
% [RT, TN, D, A, RV, RD, MD, RW] = pa_spk_readdat(fname)
% 
% Load:
%     reaction time   = RT;
%     Trialnum        = TN;
%     duration        = D;
%     Attenuation     = A;
%     ripplevelocity  = RV;
%     rippledensity   = RD;
%     modulationdepth = MD;
%     reward          = RW;
%
%

% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Initialization
fname		= pa_fcheckext(fname,'.dat');

%% BEHAVIOR
% [TN, D, A, RV, RD, MD, RT, RW] = textread (cfg.behaviorfile,...
%         '%d %d  %f %f %f %f %d %d','headerlines',3,'delimiter',';');
%     Trialnum        = TN;
%     duration        = D;
%     Attenuation     = A;
%     ripplevelocity  = RV;
%     rippledensity   = RD;
%     modulationdepth = MD;
%     reac            = RT;
%     reward          = RW;

% [TN, D, A, RV, RD, MD, RT, RW] = textread(fname,...
%         '%d %d  %f %f %f %f %d %d','headerlines',3,'delimiter',';');

fid = fopen(fname);
C	= textscan(fid,'%d %d  %f %f %f %f %d %d','headerlines',3,'delimiter',';');
[TN, D, A, RV, RD, MD, RT, RW] = deal(C{:});
fclose(fid);

	% if length(data)==(length(reac)-1)
% 	reac = reac(1:end-1);
% 	n = length(data);
% 	str = ['# Trials in SRCfile: ' num2str(n)];
% 	disp(str);
% 	n = length(reac);
% 	str = ['# Trials in DATfile: ' num2str(n)];
% 	disp(str);
% elseif length(data)==(length(reac)+1)
% 	data = data(1:end-1);
% 	n = length(data);
% 	str = ['# Trials in SRCfile: ' num2str(n)];
% 	disp(str);
% 	n = length(reac);
% 	str = ['# Trials in DATfile: ' num2str(n)];
% 	disp(str);
% elseif length(data)~=length(reac)
% 	n = length(data);
% 	str = ['# Trials in SRCfile: ' num2str(n)];
% 	disp(str);
% 	n = length(reac);
% 	str = ['# Trials in DATfile: ' num2str(n)];
% 	disp(str);
% 	error('Number of trials in src does not match number of trials in dat');
% end
