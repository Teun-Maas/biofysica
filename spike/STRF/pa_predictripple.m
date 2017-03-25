function [Corr,Y] = pa_predictripple(strfFname,dname,varargin)

% PA_PREDICTVOCAL
%
% Assumptions: 0.25 octave, 2.5 octave bandwidth

%% Initialization
close all

if nargin<1
	strfFname = 'joe6715c01b00.mat';
	% 	strfFname = 'nep';
end
if nargin<2
	dname = 'E:\DATA\Cortex\Test\Both';
end
dspFlag       = pa_keyval('display',varargin);  % sample frequency of ripple stimulus
if isempty(dspFlag)
	dspFlag = 1;
end
msgid	= 'stats:nlinfit:IllConditionedJacobian';
warning('off',msgid);

%% Reload data
cd(dname);
strfFname = pa_fcheckexist(strfFname);
load(strfFname);
[Y,Corr] = getdata(spikeP,dspFlag);


function [Y,Corr] = getdata(Spike,dspFlag)
dt = 12.5;
di = 1;
strf		= pa_spk_ripple2strf(Spike,'comp',1);
strf		= strf.strf;
%% Obtain firing rate during different periods
A	= [Spike.spiketime];
A	= A-300; % Remove the default 300 ms silent period
Dstat		= NaN(size(Spike));
Drand		= Dstat;
Velocity	= NaN(size(Spike));
Density		= Dstat;
MD			= Dstat;
for ii = 1:length(Spike)
	Dstat(ii)		= Spike(ii).stimvalues(3);
	Drand(ii)		= Spike(ii).stimvalues(4);
	Velocity(ii)	= Spike(ii).stimvalues(5);
	Density(ii)		= Spike(ii).stimvalues(6);
	MD(ii)			= Spike(ii).stimvalues(7);
end
% A			= A-Dstat;
% Some rounding, or else Matlab will not recognize certain values correctly
% due to rounding precision
Velocity	= round(Velocity*1000)/1000;
Density		= round(Density*1000)/1000;
MD			= round(MD);

% And get unique values
uV			= unique(Velocity);
uD			= unique(Density);
uM			= unique(MD);

%% Get STRF
Corr	= squeeze(NaN([numel(uV) numel(uD) numel(uM)]));
for ii = 1:numel(uV)
	for jj = 1:numel(uD)
		for kk = 1:numel(uM)
			
			selv	= Velocity==uV(ii);
			seld	= Density==uD(jj);
			selm	= MD==uM(kk);
			sel		= selv & seld & selm;
			s		= Spike(sel);
			nrep	= sum(sel);
			
			Y		= pa_genripple(uV(ii),uD(jj),uM(kk),1000,500,'display',0);
			S		= pa_spk_predict(strf,Y,dt,'display',dspFlag,'interp',di,'abs',1);
			
			tim		= [s.spiketime];
			x		= 0:(dt/di):2500;
			fr		= hist(tim,x)/(dt/di)*1000/nrep;
			
			
			sel = x>=300;
			fr = fr(sel);
				n = length(S);
				fr = fr(1:n);
			t = (0:length(fr)-1)*dt;
			r = corrcoef(S,fr);
			Corr(ii,jj,kk) = r(2);
			if dspFlag
				figure(1)
				cla;
				hold on
				plot(t,fr,'r-','LineWidth',2);
				str = ['r = ' num2str(r(2),2)];
				title(str);
				pa_verline(500,'r-');
				pa_verline(1500,'r-');
% 				ax(1) = subplot(324);
% 				ax(2) = subplot(326);
% 				linkaxes(ax,'x');
				pause;
			end
			% 			pause(1);
			% return
			% 			FR
		end
	end
end


