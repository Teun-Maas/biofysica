function data = nirs_optodedistance(data)
% DATA = NIRS_OPTODEDISTANCE(DATA)
%
% Determine the distance between fibers for every Receiver-Transmitter
% channel, and insert in DATA.opto.fiberdistance.
%
% Also corrects for a current bug in Artinis' fieldtrip software, which mislabels fibers. 
%
% See also NIRS_RCS, NIRS_OPTODELAYOUT


%% Initialization
if nargin<1
	w			= what('LR-04-2015-06-17_active');
	DataFolder	= w.path;
	cd(DataFolder)
	fname		= 'data_trans.mat';
	load(fname);
	data = data_trans;
end

%% Relevant parameters
label	= data.label; % transformed channel label, combination Receiver and Transmitter
flabel	= data.opto.fiberlabel; % fiber label
fpos	= data.opto.fiberpos; % fiber positions

label	= label(1:2:end); % deoxy labeling removed
npos	= numel(label);

xf		= fpos(:,1);
yf		= fpos(:,2);


%% Correcting a flaw in fieldtrip Artinis plugin
% for tIdx = 1:16
% 	flabel{end-16+tIdx} = ['Tx' num2str(tIdx)]; % number from 1 to 16
% end
% data.opto.fiberlabel = flabel;

%% determine distance between Receiver and Transmitter fibers
d				= NaN(npos,1);
for posIdx		= 1:npos
	str			= label{posIdx};
	chanRstr	= str(1:4);
	chanTstr	= str(end-10:end-7);
	if strcmp(chanTstr(1),'-')
		chanTstr	= chanTstr(2:end);
	end
	idxR		= strfind(flabel,chanRstr); % cell array with empty 
	idxR		= find(not(cellfun('isempty', idxR)));
	idxT		= strfind(flabel,chanTstr);
	idxT		= find(not(cellfun('isempty', idxT)));
	if strcmp(chanTstr,'Tx1') % Tx1 also appears in Tx10, Tx11, etc
		idxT	= idxT(1); % so choose the first one
	end	
	d(posIdx)	= sqrt( (xf(idxR)-xf(idxT)).^2+(yf(idxR)-yf(idxT)).^2 ); % Pythagorean theorem
end

d2				= [d d]'; d2	= d2(:); % Double for oxy and deoxy
data.opto.fiberdistance = d2; % Insert into structure

