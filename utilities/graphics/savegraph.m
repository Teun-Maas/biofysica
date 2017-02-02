function savegraph(file,type)
% SAVEGRAPH(FILE,TYPE)
%
% Functions for opening and saving graphics that operate the same for
% Windows and Macintosh operating systems.
%
% See also OPENGRAPH

%% Initialization
if nargin<1
	file = 'saveGraphOutput';
end
if nargin<2
	type = 'eps';
end

%% Maximum size
fig = gcf;
set(fig,'units','normalized','outerposition',[0 0 1 1])
fig.PaperPositionMode = 'auto';

%% Save
file	= fcheckext(file,type);
formats = {'png','jpeg','jpg','tiff','bmp','pdf','eps'};
if any(strcmp(type,formats))
	sptype = ['-d' type];
	switch type
		case 'jpg'
			sptype = '-djpeg';
		case 'bmp'
			sptype = '-dbmp16m'; %24 bit BMP file format
		case 'eps'
			sptype = '-depsc';
	end
	switch type
		case 'eps'
			if ispc && verLessThan('matlab', '8.1.0') % Windows OS for older Matlab versions
				print(sptype,'-r300','-painter',file); % to really save in Vector format
			else % Mac OS
				print(sptype,'-r300','-painters',file); % to really save in Vector format
			end
		otherwise
			print(sptype,'-r300',file);
	end
end