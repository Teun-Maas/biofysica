function [data,cfg] = pa_spk_readsrc(cfg)
% [DATA,CFG] = PA_SPK_READSRC(FNAME)
%
% Convert a stimulus-based SRC-structure in file FNAME to a trial-based
% DAT-structure. The CFG-structure contains all non-data, accessory
% information, such as filename, and analysis methods performed on the
% data.
% 
% Fields in DAT are:
%	- spikewave
%	- spiketime
%	- stimparams
%	- stimvalues
%	- timestamp
%
% Field in CFG are:
%   - filename
%   - function
%   - previous
%
% SRC-files are stored by Jann Schnupp's BrainWare for TDT systems.
%
% The trial-order, required for comparing with behavioral data, is stored in
% field trialorder. You can order the DAT-structure according to
% trial-order like this:
% >> order = [DAT.trialorder];
% >> DAT = DAT(order);
%
% So if DAT(1).trialorder = 41, this tells you that the first presented
% trial is stored in the 41st index of the DAT structure.
%
% See also PA_SPK_READDAT, PA_SPK_CHECKTRIAL

% Deciphered code from Jan Schnupp's BrainWare
% 2011 Marc van Wanrooij
% E-mail: marcvanwanrooij@neural-code.com

%% Initialization
if nargin<1 % check whether there is any input
	fname = pa_fcheckexist([],'*.src');
	cfg.filename = fname;
end
if ischar(cfg) % if the input is a string
    fname = cfg;
    clear cfg;
    cfg.filename = fname; % assume it is a filename, and not a configuration-struct
end
if ~isfield(cfg,'filename') % But if it is a struct, and does not contain a filename
	fname = pa_fcheckexist([],'*.src'); % ask for a filename
	cfg.filename = fname;
end
previous	= cfg;
data		= struct([]);
fname		= cfg.filename;
fname		= pa_fcheckext(fname,'src');
src			= readSRCfile(fname);
sets		= src.sets;

%% Clusters
[msets,nsets]	= size(sets); %#ok<*ASGLU>
[mclus,nclus]	= size(sets(1).clusters);
% [msweep,nsweep] = size(sets(1).unassignedSpikes);
n = 0;
for ii		= 1:nsets % for each set / stimulus
	cluster = sets(ii).clusters;
	for jj = 1:nclus % for each cell cluster
		sweep = cluster(jj).sweeps;
		for kk = 1:size(sweep,2) % for each repetition
			n = n+1;
% 			n	= (ii-1)*nsweep+kk; % this would have been nice, if the
% 			number of repetitions were constant
			wv						= [sweep(kk).spikes];
			if isfield(data,'spikewave') % Check whether a cluster has already been assigned to data
				[mdata,ndata] = size(data);
				if ndata>=n
					prevspikewave	= [data(n).spikewave];
					prevspiketime	= [data(n).spiketime];
					if ~isempty(wv)
						a						= [wv.shape];
						data(n).spikewave		= [prevspikewave a];
						data(n).spiketime		= [prevspiketime [wv.time]];
					end
				else
					if ~isempty(wv)
						a						= [wv.shape];
						data(n).spikewave		= a;
						data(n).spiketime		= [wv.time];
					else
						data(n).spikewave		= [];
						data(n).spiketime		= [];
					end
					
				end
			else
				if ~isempty(wv)
					a						= [wv.shape];
					data(n).spikewave		= a;
					data(n).spiketime		= [wv.time];
				else
					data(n).spikewave		= [];
					data(n).spiketime		= [];
				end
			end
		end
	end
end

%% Unassigned
n			= 0;
timestamp	= [];
for ii		= 1:size(sets,2)
	if ~isempty(sets(ii).unassignedSpikes)
		unspikes = sets(ii).unassignedSpikes;
		for jj = 1:size(unspikes,2)
			n = n+1;
% 			n	= (ii-1)*nsweep+jj;
			ts						= [unspikes(jj).timeStamp];
			timestamp				= [timestamp ts]; %#ok<AGROW>
			wv						= [unspikes(jj).spikes];
			if isfield(data,'spikewave') % Check whether a cluster has already been assigned to data
				[mdata,ndata] = size(data);
				if ndata>=n
					prevspikewave	= [data(n).spikewave];
					prevspiketime	= [data(n).spiketime];
					if ~isempty(wv)
						a						= [wv.shape];
						data(n).spikewave		= [prevspikewave a];
						data(n).spiketime		= [prevspiketime [wv.time]];
					end
				else
					if ~isempty(wv)
						a						= [wv.shape];
						data(n).spikewave		= a;
						data(n).spiketime		= [wv.time];
					else
						data(n).spikewave		= [];
						data(n).spiketime		= [];
					end
					
				end
			else
				if ~isempty(wv)
					a						= [wv.shape];
					data(n).spikewave		= a;
					data(n).spiketime		= [wv.time];
				else
					data(n).spikewave		= [];
					data(n).spiketime		= [];
				end
			end
			data(n).timestamp		= ts;
			data(n).stimparams		= sets(ii).stim.paramName;
			data(n).stimvalues		= sets(ii).stim.paramVal;
		end
	end
end
% %%  Sometimes trial/data is missing
% 
% k = [];
% for ii = 1:n
% 	if isempty(data(ii).spikewave)
% 		k = [k ii]; %#ok<AGROW>
% 	end
% end
% indx	= 1:size(data,2);
% sel		= ~ismember(indx,k);
% indx	= indx(sel);
% data	= data(indx);
% n		= size(data,2);
% whos data

%% Aborted trials
timestamp	= [data.timestamp];
mn			= min(timestamp);
timestamp	= (timestamp-mn);
abort = 0;
if isfield(src,'comments');
	if isfield(src.comments,'timeStamp');
		abort = 1;
		aborted_ts		= [src.comments.timeStamp]; % timestamps for abortions (not the same timestamps as for trials)
		aborted_real	= char({src.comments.sender}); % Sender = ABORT SIGNAL if trial is aborted
		indx			= strmatch('ABORT SIGNAL',aborted_real);
		aborted_ts		= aborted_ts(indx);
		aborted_ts		= aborted_ts-mn;
		trial_abort		= aborted_ts;
		% Somehow ABORTS line up with RW if we do the following
		for ii = 1:length(aborted_ts)
			dff				= aborted_ts(ii)-timestamp;
			[~,indx]		= min(abs(dff)); % minimal absolute difference
			trial_abort(ii) = indx+1; % and add 1 trial
		end
	end
end

%% Trial Order in experiment
[~,trialorder]	= sort(timestamp);
for ii = 1:n
	data(ii).trialorder = trialorder(ii);
	b		= data(ii).spiketime;
	sz		= size(b);
	data(ii).trial = repmat(ii,sz);
	if abort
		if ismember(ii,trial_abort)
			data(ii).aborted = 1;
		else
			data(ii).aborted = 0;
		end
	else
		data(ii).aborted = 0;
	end
end

% %% Sort by trial
% data = data([data.trialorder]);

%% Bookkeeping
% Add function name to configuration
cfg.function = mfilename('fullpath');

% Add input/previous configuration to cfg structure
cfg.previous = previous;

function data=readSRCfile(fname)
%********************************
% function data=readSRCfile(fname);
%  reads Spike ReCording file "fname" generated with BrainWare.
%  and generates a structured array "data" which mimicks the 
%  data object hierarchy used by BrainWare to organise spike data internally
%
data=[];
f=fopen(fname,'r');
while ~feof(f)
    data=[data readDataObject(f)]; %#ok<AGROW>
end;
fclose(f);

function dataObject=readDataObject(f)
%************************************
objId=fread(f,1,'uint16'); % get ID number for next data object
if ~isempty(objId) % i.e. while file f has not run to end of file
    switch objId
    case 29079, % Fixed length (40 time bin) action potential
        dataObject=readFixedLenSpike(f);
    case 29115, % Variable length action potential
        dataObject=readVarLenSpike(f);
    case 29113; % display information object, 34 bytes wide - skip
        dataObject=[]; fread(f,34,'uint8'); 
    case 29100; % old style display information object, 4 bytes wide - skip
        dataObject=[]; fread(f,4,'uint8'); 
    case 29093; % old style data set collection, same as list
        dataObject=readList(f);
    case 29112; % data set collection with comments
        dataObject=readSetCollectionObject(f);
    case 29114; % variable ADperiod data set collection 
        dataObject=readVarPeriodSetCollectionObject(f);
    case 29117; % V8 data set collection 
        dataObject=readV8SetCollectionObject(f);
    case 29110; % time stamped data sweep 
        dataObject=readTimeStampedSweep(f);
    case 29082; % normal data sweep 
        dataObject=readSweep(f);
    case 29106; % individual data set
        dataObject=readDatasetObj(f);
    case 29083; % collection (list) of data objects
        dataObject=readList(f);
    case 29084; % spike record (i.e. spike data for a cluster, but without cluster bounds
        dataObject=readSpkRec(f);
    case 29107; % old style cluster with non-standard 6 byte float boundaries
        dataObject=readClusterV1(f);
    case 29116; % new style cluster with IEEE 4 byte float boundaries
        dataObject=readClusterV2(f);
    case 29091; % collection of clusters - same structure as list
        dataObject=readList(f);
    case 29109; % stimulus descriptor object
        dataObject=readStimObj(f);
    otherwise 
        fclose(f);
        error(['read unrecognized data object ID ' int2str(objId)]);
    end;
else 
    dataObject=[];
end; 
% done function readDataObject


function data=readVarPeriodSetCollectionObject(f)
%*************************************************
anADperiod=fread(f,1,'float32'); % DA conversion clock period in microsec
data=readSetCollectionObject(f);
data.ADperiod=anADperiod;
% done function readVarPeriodSetCollectionObject


function data=readV8SetCollectionObject(f)
%*************************************************
data=readVarPeriodSetCollectionObject(f);
data.sortInfo=readSorter(f);
% dummy=fread(f,1,'int16'); % Commented by Marc


% done function readV8SetCollectionObject


function data=readSorter(f)
%*************************************************
data.nTimeSlices=fread(f,1,'int16');
data.timeslice=[];
for ii=1:data.nTimeSlices,
    slice=[];
    fread(f,1,'int16');
    slice.maxValid=fread(f,1,'double');
    slice.nClust=fread(f,1,'int16');
    slice.cluster=[];
    for jj=1:slice.nClust,
        clust=[];
        dummy2=fread(f,1,'int16'); %#ok<NASGU> % variable not used, but needs to be read
        clust.numChans=fread(f,1,'int16');
        clust.elliptic=fread(f,10*clust.numChans,'uint8'); % skip boolean values indicating elliptic feature boundary dimensions
        clust.boundaries=fread(f,20*clust.numChans,'float32'); % read boundaries
        slice.cluster=[slice.cluster, clust];
    end;
    data.timeslice=[data.timeslice, slice];
end;
    
% done function readV8SetCollectionObject




function data=readSetCollectionObject(f)
%****************************************
data.NChannels=fread(f,1,'uint8'); % read number of electrode channels in set
data.sets=readList(f); % read individual datasets
% data.side=char(fread(f,1,'char')); % read "side of brain" info
data.side=fread(f,1,'*char'); % read "side of brain" info % changed by Marc

% read comments
numComments=fread(f,1,'int16');
data.comments=[];
for ii=1:numComments,
    data.comments=[data.comments readCommentObj(f)];
end;
% done function readSetCollectionObject


function aList=readList(f)
%**************************
numElements=fread(f,3,'int16'); %read three integers. 
% The 1st integer gives number of data sets. 
% (The others can be ignored).
aList=[];
for ii=1:numElements(1)
    % each element is read in turn by readDataObject() and appended to the list
    newElement=readDataObject(f);
    aList=[aList newElement]; %#ok<AGROW>
end;

function set=readDatasetObj(f)
%****************************** 
% read and return the dataset object
set.stim=readDataObject(f); % read the stimulus for this set
set.unassignedSpikes=readDataObject(f); % read sweep collection of unassigned spikes
set.clusters=readDataObject(f); % read collection of spike clusters
set.sweepLen=fread(f,1,'int32'); % read sweep length for dataset
% done readDatasetObj()


function stim=readStimObj(f)
%****************************
% read a stimulus object from file
stim.numParams=fread(f,1,'int16'); % read number of stimulus parameters
stim.paramName=[]; % initialise list of parameter names
% read stimulus parameter names in turn
for ii=1:stim.numParams,
    fread(f,1,'char'); % skip one byte
    paramNameLen=fread(f,1,'uint8'); % read length of next parameter name
    nextName=fread(f,paramNameLen,'char')'; % read next parameter name 
    nextName=cellstr(char(nextName));% ... convert to string
    stim.paramName=[stim.paramName, nextName]; % ...and add it to the list of parameter names
end;
% finally read parameter values - an array of numParams 32-bit floating point numbers
stim.paramVal=fread(f,stim.numParams,'float32');
% done readStimObj()


function sweep=readTimeStampedSweep(f)
%*************************************
sweep.timeStamp=fread(f,1,'float64'); % read timeStamp (number of days since dec 30th 1899)
sweep.spikes=readList(f); % read spike list
% done readTimeStampedSweep(f);

function sweep=readSweep(f)
%***************************
sweep.spikes=readList(f); % read spike list
% done readSweep(f);

% readClusters  does not seem to be used, so removed by Marc
% function clusters=readClusters(f)
% %*********************************
% objId=fread(f,1,'uint16');
% if objId ~= 29091,
%     error(['Trying to read data cluster collection but found bad object id: ' int2str(objId)]);
% end;
% numClust=fread(f,3,'int16'); % read three integers. 1st int is number of clusters.
% clusters=[];  % initialise list of clusters
% for ii=1:numClust;
%     clusters=[clusters readCluster(f)]; % read individual clusters and add to list
% end;
% % done readClusters(f);


function cluster=readSpkRec(f)
%********************************
fread(f,1,'uint16'); % skip two bytes
% read ID string (cluster 'name')
numChars=fread(f,1,'uint16');
anIDStr=fread(f,numChars,'char')';
cluster.IdString=char(anIDStr);
cluster.sweepLen=fread(f,1,'int32'); % sweep length in ms
cluster.respWin=fread(f,4,'int32'); % response and spon period boundaries
cluster.sweeps=readDataObject(f); % data sweeps
% done readSpkRec(f);


function cluster=readClusterV1(f)
%********************************
cluster=readSpkRec(f);
% cluster boundaries are stored in 48-bit floating point numbers which is not supported in Matlab - skip
fread(f,18*6,'uint8'); % skip boundaries
cluster.elliptic=fread(f,9,'uint8'); % skip boolean values indicating elliptic feature boundary dimensions
% done readClusterV1(f);



function cluster=readClusterV2(f)
%********************************
fread(f,1,'uint16'); % skip two bytes
% read ID string (cluster 'name')
numChars=fread(f,1,'uint16');
anIDStr=fread(f,numChars,'char')';
cluster.IdString=char(anIDStr);
cluster.sweepLen=fread(f,1,'int32'); % sweep length in ms
cluster.respWin=fread(f,4,'int32'); % response and spon period boundaries
cluster.sweeps=readDataObject(f); % data sweeps
% cluster boundaries stored in IEEE 32-bit floats
cluster.boundaries=fread(f,18,'float32'); % read boundaries
cluster.elliptic=fread(f,9,'uint8'); % skip boolean values indicating elliptic feature boundary dimensions
% done readClusterV2(f);


function spike=readFixedLenSpike(f)
%***********************************
spike.time=fread(f,1,'float32'); % read spike time stamp in ms since start of sweep
spike.shape=fread(f,40,'int8'); % read spike shape;
spike.trig2=fread(f,1,'uint8'); % point of return to noise
% done readFixedLenSpike


function spike=readVarLenSpike(f)
%**********************************
numPts=fread(f,1,'uint8');
spike.time=fread(f,1,'float32'); % read spike time stamp in ms since start of sweep
spike.shape=fread(f,numPts,'int8'); % read spike shape;
spike.trig2=fread(f,1,'uint8'); % point of return to noise
% done readVarLenSpike


function comment=readCommentObj(f)
%**********************************
comment.timeStamp=fread(f,1,'float64'); % read timeStamp (number of days since dec 30th 1899)
numChars=fread(f,1,'int16');
comment.sender=char(fread(f,numChars,'char')');
numChars=fread(f,1,'int16');
comment.text=char(fread(f,numChars,'char')');
% done readCommentObj()
