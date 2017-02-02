function dbhl = pa_dba2dbhl(dba,f,method)
% DBHL = PA_DBA2DBHL(DBA)
%
% Convert dBA to dBHL
%
% As determined from figure 6 in:
% 1) http://www.batod.org.uk/content/resources/audiology/refreshers/general/G3-decibels.pdf
% Also:
% 2) http://ec.europa.eu/health/opinions/en/hearing-loss-personal-music-player-mp3/figtableboxes/table-2.htm
% 3) http://thebsa.org.uk/docs/Guidelines/soundfieldguidelinesfeb2008.pdf
%
% See also PA_AWEIGHT, PA_TDT_BOSE_AUDIOGRAM

% 25 January 2013 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

if nargin<2
	f = 100:100:10000;
end
if nargin<1
	dba = zeros(size(f));
end
if nargin<3
	method = 1;
end

switch method
	case 1
		% 1) http://www.batod.org.uk/content/resources/audiology/refreshers/general/G3-decibels.pdf
		freq = [250 500 1000 2000 4000];
		a2hl = [-20 -17 -10 -10 -10];
	case 2
		% 2) http://ec.europa.eu/health/opinions/en/hearing-loss-personal-music-player-mp3/figtableboxes/table-2.htm
		freq = [250 500 1000 2000 4000 8000];
		a2hl = -[12 5 2 -2 -5 13];
	case 3
		% 3) http://thebsa.org.uk/docs/Guidelines/soundfieldguidelinesfeb2008.pdf
		% for this method we need to convert dBA to dB SPL first
		freq = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
		W = pa_aweight(f);
		dba = dba-W;
		a2hl = -[2.4 0.3 -0.1 0.5 -0.1 -0.3 -2.8 -4 -0.5 4.4];
	case 4
		%% Average
		freq = [250 500 1000 2000 4000];
		a2hl = [-20 -17 -10 -10 -10];
		ai		= interp1(freq,a2hl,f,'cubic','extrap');
		dbhl1	= dba+ai;
		
		freq = [250 500 1000 2000 4000 8000];
		a2hl = -[12 5 2 -2 -5 13];
		ai		= interp1(freq,a2hl,f,'cubic','extrap');
		dbhl2	= dba+ai;
		
		% for this method we need to convert dBA to dB SPL first
		freq = [250 500 750 1000 1500 2000 3000 4000 6000 8000];
		W = pa_aweight(f);
		dba = dba-W;
		a2hl = -[2.4 0.3 -0.1 0.5 -0.1 -0.3 -2.8 -4 -0.5 4.4];
		ai		= interp1(freq,a2hl,f,'cubic','extrap');
		dbhl3	= dba+ai;
		dbhl = (dbhl1+dbhl2+dbhl3)/3;
		return
end

ai		= interp1(freq,a2hl,f,'cubic','extrap');
dbhl	= dba+ai;

%%
% plot(freq,a2hl,'k.-')
% hold on
% plot(f,ai,'ro-')
