function [data,x,y,label] = nirs_optodelayout(fname,varargin)
% [X,Y,LABEL] = NIRS_OPTODELAYOUT(FNAME)
%
% Get optode layout [X,Y] for deep channels LABEL for Kennan experiment
% BEFORE NIRS_RCS
%
% FNAME is a mat-filename obtained with fieldtrip NIRS Artinis plugin, or a
% data stuct
%
% See also NIRS_RCS

%% Initialization
if nargin<1
	w			= what('LR-04-2015-06-17_active');
	DataFolder	= w.path;
	cd(DataFolder)
	fname		= 'data.mat';
end
if ischar(fname)
	load(fname)
elseif isstruct(fname)
	data	= fname;
end
dispFlag		= keyval('disp',varargin,false);
layFlag			= keyval('lay',varargin,false);
% f				= keyval('fig',varargin,2);


%% Optode distance
if ~isfield(data.opto,'fiberdistance')
	data	= nirs_optodedistance(data);
end
% data = nirs_optodedistance(data);


%% Relevant parameters
pos		= data.opto.chanpos; % channel position
label	= data.label; % transformed channel label
flabel	= data.opto.fiberlabel; % fiber label
fpos	= data.opto.fiberpos; % fiber positions
d		= data.opto.fiberdistance; % fiber distance
x		= pos(:,1);
y		= pos(:,2);

%% set anterior to top
[x,y,~] = yaw(x,y,ones(size(x)),90);

%% mirror right hemisphere (so that dorsal is inside, and ventral outside)
x		= x+0.5;
sel		= x>0;
x		= x-0.5;
x(sel)	= -x(sel);
x(sel)	= x(sel)-min(x(sel));

%% Same matrix manipulations for fibers (instead of channels)
xf		= fpos(:,1);
yf		= fpos(:,2);
nfpos	= numel(xf);

[xf,yf,~]	= yaw(xf,yf,ones(size(xf)),90);
xf			= xf+0.5;
sel			= xf>0;
xf			= xf-0.5;
xf(sel)		= -xf(sel);
xf(sel)		= xf(sel)-min(xf(sel));

mux			= mean(xf);
muy			= mean(yf);
xf			= xf-mux;
yf			= yf-muy;
x			= x-mux;
y			= y-muy;

%% Separate left and right hemisphere
sel			= xf<0;
xf(sel)		= xf(sel)-6;
xf(~sel)	= xf(~sel)+6;

sel			= x<0;
x(sel)		= x(sel)-6;
x(~sel)		= x(~sel)+6;

%% Correcting a flaw in nirs
% No longer necessary? Sept 30 2016
% for tIdx = 1:16
% 	flabel{end-16+tIdx} = ['Tx' num2str(tIdx)];
% end

data.chanpos = [x y];


%% lay.mat
% pos		= data.opto.chanpos; % channel position
% label	= data.label; % transformed channel label
% flabel	= data.opto.fiberlabel; % fiber label
% fpos	= data.opto.fiberpos; % fiber positions
if layFlag
	
	layd		= data.opto.fiberdistance; % fiber positions
	
	sel			= layd==max(layd);
	lay.pos		= [x(sel) y(sel)]*50;
	lay.width	= repmat(100,size(lay.pos,1),1);
	lay.height	= repmat(100,size(lay.pos,1),1);
	lay.label	= data.label(sel);
	cd('/Users/marcw/Gitlab/biofysica/nirs/layout');
	save('ciavlay.mat','lay');
end
%%

%% Graphics
if dispFlag(1)
	
	x		= x(1:2:end)
	y		= y(1:2:end)
	label	= label(1:2:end);
	npos	= numel(x);
	keyboard
	d		= d(1:2:end);
	% Deep vs shallow channel
	ud		= unique(d);
	shallow = min(ud);
	deep	= max(ud);
	
	hold on
	% plot(x,y,'ko','MarkerFaceColor','w')
	
	%%
	for posIdx = 1:npos
		posIdx
		if d(posIdx)==shallow
			col = 'b';
		elseif d(posIdx)==deep
			col = 'r';
		else
			col = 'k';
		end
		
		str = label{posIdx};
		str = str(1:end-7);
		text(x(posIdx),y(posIdx),str,'HorizontalAlignment','center','Color',col)
	end
	
	%%
	for posIdx = 1:nfpos
		str = flabel{posIdx};
		text(xf(posIdx),yf(posIdx),str,'HorizontalAlignment','center','Color',[.7 .7 .7])
	end
	axis([-25 25 -25 25]);
	% axis square;
	xlabel('left < X (cm) > right');
	ylabel('posterior < Y (cm) > anterior');
	text(0,0,'dorsal','horizontalalignment','center','verticalalignment','middle','rotation',90);
	
	set(gca,'XTick',-24:1.5:24,'YTick',-24:1.5:24);
	grid on
end

%% .lay
% The layout file is a plain ascii file with the extention *.lay
% The 1st column is the channel number in the layout file, it is not used any more by the plotting functions, but should be present in the layout file.
% The 2nd and 3rd are the X-position and Y-position.
% The 4th and 5th column specify the width and height The 6th column is a string with the channel label.
% if layFlag
% 	fid = fopen('kennanhelmet.lay','w');
%
% 	for ii = 1:numel(x)
% 		fprintf(fid,'%i\t%i\t%i\t%i\t%i\t%s\n',ii,x(ii),y(ii),w(ii),h(ii),label{ii});
% 	end
% 	for ii = 1:numel(x)
% 		fprintf(fid,'%i\t%i\t%i\t%i\t%i\t%s\n',ii,x(ii),y(ii),w(ii),h(ii),delabel{ii});
% 	end
% 	fclose(fid);
% end
