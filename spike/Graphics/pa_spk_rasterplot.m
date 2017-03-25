function h = pa_spk_rasterplot(T,varargin)
% PA_SPK_RASTERPLOT(T)
%
% Display raster plot of spikes in structure T. T should contain spike
% timings for every trial in the field spiketime.
%
% Alternatively, T can be a MxN length vector, with M t
%at times T (in samples) for NTRIAL trials,
% each of length TRIALLENGTH samples, sampling rate = 1kHz. SpikeT are
% hashed by the trial length. 
%
% PA_SPK_RASTERPLOT(T,NTRIAL,TRIALLENGTH,FS)
% Plots the rasters using sampling rate of FS (Hz)
%
% H = PA_SPK_RASTERPLOT(T,NTRIAL,TRIALLENGTH,FS)
% Get handle of rasterplot figure
%
% PA_GETPOWER(...,'PARAM1',val1,'PARAM2',val2) specifies optional
% name/value pairs. Parameters are:
%	'color'	- specify colour of graph. Colour choices are the same as for
%	PLOT (default: k - black).
%
%  Example:
%		Ntrials = 50;
%		Ltrial  = 1000;
%		nspikes = 1000;
%		T       = pa_spk_genspiketimes(nspikes,Ntrials*Ltrial);
%       h       = pa_spk_rasterplot(T,'Ntrials',Ntrials,'Ltrial',Ltrial);
%
% More information: 
%		"Spikes - exploring the neural code", Rieke et al. 1999, figure	2.1
%		"Matlab for Neuroscientists", Wallisch et al. 2009,  section 13.3.1

% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Initialization
Ntrials         = pa_keyval('Ntrials',varargin);
if isempty(Ntrials)
	Ntrials		= 1;
end
Ltrial         = pa_keyval('Ltrial',varargin);
if isempty(Ltrial)
	Ltrial		= round(length(T)/Ntrials);
end
Fs         = pa_keyval('Fs',varargin);
if isempty(Fs)
	Fs			= 1000; % (Hz)
end

%% Plot variables
col         = pa_keyval('color',varargin);
if isempty(col)
	col			= 'k';
end

%% Plotting
if isstruct(T) % if we have a structure containing spike timings
	Ntrials = length(T);
	hold on;
	for ii = 1:Ntrials
		t = T(ii).spiketime;
		t = t/Fs*1000; % spike times in (ms)
		for jj = 1:length(t)
			line([t(jj) t(jj)],[ii-1 ii],'Color',col);
		end
	end
	% Ticks and labels
	ylim([0 Ntrials]);
	set(gca,'TickDir','out');
	xlabel('Time (ms)');
	ylabel('Trial #');
	% It is inconvenient to get a line handle for EVERY spike, therefore we just
	% give the handle of the current axis
	h = gca;
else % if we have just on- and offsets
	% plot spikes
	trials		= ceil(T/Ltrial);	% get trial number for each spike time
	time		= mod(T,Ltrial);	% get time in trial
	time(~time) = Ltrial;			% Time 0 is actually last sample in trial
	numspikes	= length(T);
	x			= NaN(3*numspikes,1);
	y			= NaN(3*numspikes,1);
	
	% Trials
	y(1:3:3*numspikes)		= (trials-1);
	y(2:3:3*numspikes)		= y(1:3:3*numspikes)+1;
	
	% Time scale
	x(1:3:3*numspikes)		= time*1000/Fs; % in ms
	x(2:3:3*numspikes)		= time*1000/Fs;
	
	h = plot(x,y,'Color',col);
	xlim([1 Ltrial*1000/Fs]);
	ylim([0 Ntrials]);
	
	% Ticks and labels
	set(gca,'TickDir','out');
	xlabel('Time re Onset (ms)');
	ylabel('Trial #');
end
