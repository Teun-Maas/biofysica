function [nirs_data, events] = oxysoft2matlab(filename, target, savename, verbose, projevents)
% This function reads in data exported from the Artinis software Oxysoft
% and returns a matlab-structure in a format that can be imported by other
% Matlab toolboxes, such as NIRS-SPM or Homer. The export can also be
% saved if desired.
%
% Use as:
%    nirs_data = oxysoft2matlab(filename, target, savename, verbose)
% where
%   filename -   string or cell, location of the file to be imported. The
%                 file has to originate from Oxysoft and must be an
%                 .oxy3-file, or an .oxyproj-file. If specifying an .oxyproj
%                 file, all .oxy3-files stored in there are imported.%
%                 If using a cell-array and mixing oxy3 and oxyproj files,
%                 only the common oxy3-files are imported.
%   target   -   string, the target toolbox. Can be 'rawOD', 'oxy/dxy',
%                 'nirs-spm', 'homer', or 'nap'.
%   savename -   string or cell, location where the export should be saved.
%                 Please specify without extension - it is automatically
%                 added. If empty, data will not be saved. If a cell, has
%                 to be of the same length as filename, which also has to
%                 be a cell-array then.
%   verbose  -   boolean, true (default) means text output in command
%                 window.
%
% All input arguments are optional. UI dialogs will open and ask for the
% necessary information (what data to convert to which toolbox and where to
% save) when an input argument is missing.
% Data will be stored in the format of the toolbox, such that you can open
% the data from within the toolbox. The output argument of this function
% contains the data in the toolbox-specific format. If you specify 'raw'
% export format, a data structure with optical densities and positions of
% the system will be returned. You can transform the data to oxy- and
% deoxyvalues using the function transformoxy3od.m.
%
% An export of 'rawOD' will be returned as:
%         nirs_data.OD           - matrix NxM, where N=#samples and
%                                   M=#measurements (transmitter wavelengths
%                                   and receiver combinations). Contains the
%                                   raw optical densities.
% An export of 'oxy/dxy' values will be returned as:
%         nirs_data.oxyvals      - matrix NxM, where N=#samples and
%                                   M=#measurements (transmitter wavelengths
%                                   and receiver combinations). Contains the
%                                   relative change in oxygenated blood.
%         nirs_data.dxyvals      - matrix NxM, where N=#samples and
%                                   M=#measurements (transmitter wavelengths
%                                   and receiver combinations). Contains the
%                                   relative change in deoxygenated blood.
% The export to 'oxy/dxy' can also contain additional fields for TSI data.
% These are absO2Hb, absHHb (absolute concentrations of (de-)oxygenated
% blood), TSI (tissue saturation index) and TSI_FF (TSI fit factor).
%
% Additional fields in the data structure for an export to rawOD or oxy/dxy
% data are:
%         nirs_data.Fs           - scalar, sampling frequency of the
%                                   recording
%         nirs_data.time         - vector Nx1, , where N = #samples.
%                                   General time axis in seconds, starting
%                                   at 0s.%
%         nirs_data.DPF          - vector, differential path length per
%                                   optode template
%         nirs_data.wavelengths  - matrix LxT, where L=#wavelengths and
%                                   T=#transmitters. Contains wavelengths
%                                   per transmitter.
%         nirs_data.label        - cell 1xC, where C=#channels. Unique
%                                   channel names.
%         nirs_data.RxLabel      - vector 1xR, where R=#receivers. Unique
%                                   receiver labels.
%         nirs_data.TxLabel      - vector 1xT, where T=#receivers. Unique
%                                   transmitter labels.
%         nirs_data.distance     - vector 1xC, where C=#channels. Distance
%                                   in cm between receiver and transmitter.
%         nirs_data.chanPos      - vector 2xC, where C=#channels. Position
%                                   of the channels in cm relative in an
%                                   arbitrary coordinate system.
%         nirs_data.transPos     - vector 2xT, where T=#transmitters.
%                                   Position of the transmitters in cm
%                                   relative in an arbitrary coordinate
%                                   system.
%         nirs_data.receiPos     - vector 2xC, where R=#receivers. Position
%                                   of the receivers in cm relative in an
%                                   arbitrary coordinate system.
%         nirs_data.Rx_TxId      - vector, 2xM, where M=#measurements.
%                                   Indexes which receivers (first row) and
%                                   transmitters (second row) were used to
%                                   measure the optical densities.
%         nirs_data.ADvalues     - matrix NxA, where N=#samples
%                                   A=#additonal channels. All data of the
%                                   additional channels that were attached
%                                   to the system.
%         nirs_data.metaInfo      - struct. Meta and header information of
%                                   the data.
%
% In case of an export to only ODs ('rawOD'), the structure will additionally
% have these fields:
%         nirs_data.ODlabel      - cell 1xC, where M=#measurements.
%                                   Respective channel names and wavelength
%                                   combination. Redundant with
%                                   .wavelength, .label and .Rx_TxId
%         nirs_data.ODchanPos    - vector 2xM, where M=#measurement. Position
%                                   of the measurement in cm relative in an
%                                   arbitrary coordinate system. Redundant
%                                   with .chanPos and .Rx_TxId.
%
% In case of an export to only ODs ('rawOD'), you should be able to use the
% transformoxy3od function in the private/ directory to convert to
% concentrations. use as
%       transformoxy3od(nirs_data.OD, nirs_data.metaInfo, nirs_data.DPF)
% For more, see help of that function
%
%
% Version 1.21, copyright (c) by Artinis Medical Systems http://www.artinis.com
% Author J�rn M. Horschig, jorn@artinis.com
%

%%
% Version history
%
% v 1.21    fix - corrected OD extraction when not all Tx but more than 
% 08/24/16        one Rx was used
%
% v 1.20    fix - single oxyproj selection handled correctly again
% 08/08/16
%
% v 1.19    fix - Condition names properly exported to Homer2
% 07/18/16  fix - when optodetemplate not found in xml, ask for alternative
%           enh - consecutive numbers in multi-system optodetemplates supported (?)
%
% v 1.18    enh - Condition names properly exported to Homer2
% 07/05/16 
%
% v 1.17    fix - filename selection in UI dialog fixed
% 07/04/16  fix - mixup of triggers/events order removed
%
% v 1.16    fix - removed ~ references to increase backwards compatbility
% 06/20/16
%
% v 1.15    enh - multiple TSI system support
% 06/01/16
%
% v 1.14    fix - NAP single event bug
% 05/25/16  fix - projevents initialized
%
% v 1.13    fix - conditions when to ask for saving export
% 05/24/16  fix - project events return a double instead of a string
%           fix - oxyproj read in
%           fix - events export to Homer2
%           enh - enabling 3D digitized position read-in
%           enh - allowing NAP-export for all templates%
%           doc - timestamp fix of v 1.12
%
% v 1.12    list of improvements and bugfixes:
% 05/20/16   fix - automatically adding oxy3 extension
%            enh - allowing for batch processing of oxy3-files
%            enh - enabling reading in of oxyproj files
%            enh - allowing filename and savename to be a cell-array
%
% v 1.11    fixed a critical optodetemplate naming bug
% 05/03/16
%
% v 1.10    fixed a critical optodetemplate naming bug
% 04/25/16
%
% v 1.09    verbose documentation added, OD vals checked again
% 02/22/16
%
% v 1.08    multi-channel, raw OD channel labels corrects
% 04/08/15   timestamp of v1.07 corrected
%
% v 1.07    added information on used templates and optodetemplates.xml
% 06/07/15   also changed precedence of optodetemplates.xml
%            and bugfixes (TSI+noTSI, wavelengths of diff. subtemplates)
%
% v 1.06    added SpatialUnit ('cm') to better support Homer2 v2.0
% 12/06/15   also fixed a bug for auxiliary channels (AD channels)
%
% v 1.05    improved memory handling and increased backwards compatibitliy
% 17/02/15   and also fixed channel label for raw transformation and structure for Homer
%
% v 1.04    fixed a bug pertaining to reading out of event names and in
% 12/22/14   case of having recorded multichannel TSI
%
% v 1.03    fixed a bug and inconsistency regarding reading out of event names
% 12/08/14
%
% v 1.02    fixed another crash during parsing of the header when Version
% 11/11/14   was not present in metaInfo and extended rawOD output with
%            the (redundant) fields ODlabel and ODchanPos.
%
% v 1.01    removed the use of strjoin for backwards compatibility and fixed
% 11/04/14   a crash during parsing of the header when Version was not
%            present in metaInfo
%
% v 1.0     initial release
% 10/30/14
%
% See also test_oxysoft2matlab at Artinis private server
% Known todos:
%   - add support for .txt, .xml and a .oxy/.evt
%   - add support for 'fosa' and 'p3' (aka 'potato')
%
% other toolboxes which might be considered later:
%   NBT toolbox: www.nbtwiki.net
%
%   LIPSIA
%
%   TOAST++
%       http://web4.cs.ucl.ac.uk/research/vis/toast/
%       M. Schweiger and S. R. Arridge, The Toast++ software suite for for-
%       ward and inverse modeling in optical tomography, J Biomed Opt (2014)
%
%   nirsLAB or NAVI
%       http://www.nirx.net/software/nirslab---nirs-topography-software-by-nirx
%       http://www.nirx.net/software/navi-nirs-tomography-software-by-nirx
%
%   EasyTopo
%       https://sites.google.com/site/fenghuatian/software/easytopo

%% input error handling
if nargin < 5
  projevents = [];
  if nargin < 4
    verbose = true;
    if nargin < 3
      savename = '';
    end
  end
end

% check if all input arguments were specified, otherwise use a UI dialog
if nargin < 1 || isempty(filename)
  w = true;
  while w
    [filename,pathname] = uigetfile({'*.oxyproj; *.oxy3', 'Oxysoft Files (*.oxy3, *.oxyproj)'}, ...
      'Select one or several .oxy3/.oxyproj-file(s)', ...
      'MultiSelect', 'on');
    if isscalar(filename) && filename ==0
      if verbose
        fprintf('[Artinis] user aborted.\n');
      end
      return;
    end
    if ischar(filename)
      filename = {filename};
    end
    for f=filename
      if ~strcmp(f{1}(end-4:end), '.oxy3') && ~strcmp(f{1}(end-7:end), '.oxyproj')
        errorbox('Please only select oxy3- or oxyproj-files.')
        w = true;
        break;
      else
        w = false;
      end
    end
  end
  for f=1:numel(filename)
    filename{f} = fullfile(pathname, filename{f});
  end
end

% check whether target toolbox is supported
% supported_targets = {'nirs-spm', 'homer', 'nap', 'fosa', 'p3', 'potato'};
supported_targets = {'nirs-spm', 'homer', 'nap', 'rawOD', 'oxy/dxy'};
if nargin < 2 || isempty(target )
  target_idx = menu('Choose a target toolbox',supported_targets);
  target     = supported_targets{target_idx};
  drawnow;
  fprintf('[Artinis] Selected target toolbox %s.\n', target);
end

if ~ismember(target, supported_targets)
  tgts = sprintf('''%s'', ', supported_targets{1:end-1});
  tgts = sprintf('%s and ''%s''', tgts(1:end-2), supported_targets{end});
  error('[Artinis] target toolbox cannot be identified. please choose among %s.', tgts);
end

supported_extensions = {'oxy3', 'oxyproj'};
% supported_extensions = {'.txt', '.xml', '.oxy3', '.oxy', '.evt'};

measurements = struct('name', '', 'file', '', 'events', struct('value', -1, 'onset', inf, 'useSeconds', false), 'positions', struct('fid', struct('name', '', 'pos', []), 'rx', struct('name', '', 'pos', []), 'tx', struct('name', '', 'pos', [])), 'isoxy3', false);
if iscell(filename) && numel(filename)>1
  % check whether each file has an extension, does exist and is supported
  [filepath, filename, extension] = cellfun(@fileparts, filename, 'UniformOutput', false);
  for f=1:numel(filename)
    
    % assuming oxy3 extension if not provided
    if isempty(extension{f})
      if verbose
        fprintf('[Artinis] No file extension provided, assuming .oxy3\n');
      end
      %filename = [filename '.oxy3'];
      extension{f} = '.oxy3';
    end
    
    % check whether file exists
    if ~exist(fullfile(filepath{f}, [filename{f} extension{f}]), 'file')
      error('[Artinis] file ''%s'' cannot be located. please specify a valid filename.', [filepath{f} extension{f}]);
    end
    
    % create a 'measurement'-list of files
    switch extension{f}
      case '.oxyproj'
        % parse the content here
        fprintf('[Artinis] Reading in oxyproj-file... please wait...')
        measurements = [measurements parseoxyproj(fullfile(filepath{f}, [filename{f} extension{f}]))];
        fprintf('done!\n')
      case '.oxy3'
        % handling the * wildcard here
        if ismember('*', filename{f})
          files = dir(fullfile(filepath{f}, [filename{f} extension{f}]));
          for fl=files'
            measurements(end+1).file = [fl.name extension{f}];
            measurements(end).events = [];
            measurements(end).isoxy3 = true;
          end
        else
          measurements(end+1).file = fullfile(filepath{f}, [filename{f} extension{f}]);
          measurements(end).events = [];
          measurements(end).isoxy3 = true;
        end
        
      otherwise
        exts = sprintf('.%s-, ', supported_extensions{1:end-1});
        exts = sprintf('%s or .%s-', exts(1:end-3), supported_extensions{end});
        error('[Artinis] file extension not supported. please select a file in %sformat.', exts);
    end
  end
  
  % check for each entry int he measurement-list if an 'isoxy3' does exist
  m=1;
  while m<=numel(measurements)
    if ~measurements(m).isoxy3
      % check if another measurement with the same name exists
      midx = ismember(cellfun(@lower, {measurements(:).file}, 'UniformOutput', false), lower(measurements(m).file)) & [measurements(:).isoxy3];
      midx = find(midx);
      if isempty(midx)
        % delete this measuement
        measurements(m) = [];
      else
        % join the event of this measurement with the other one
        for mi = midx
          measurements(mi).events = [measurements(mi).events measurements(m).events];
          measurements(mi).positions = [measurements(mi).positions measurements(m).positions];
        end
        m = m+1;
      end
    else
      m = m+1;
    end
  end
  
  % and remove all those which are of non-oxy3 origin
  midx = [measurements.isoxy3]==false;
  measurements(midx) = [];
else
  if iscell(filename)
    filename = filename{1};
  end
  % handle a single file here
  [filepath, filename, extension] = fileparts(filename);
  
  % assuming oxy3 extension if not provided
  if isempty(extension)
    if verbose
      fprintf('[Artinis] No file extension provided, assuming .oxy3\n');
    end
    %filename = [filename '.oxy3'];
    extension = '.oxy3';
  end
  
  % check whether file exists
  if ~exist(fullfile(filepath, [filename extension]), 'file') && ~ismember('*', filename)
    error('[Artinis] file ''%s'' cannot be located. Please specify a valid filename.', fullfile(filepath, [filename extension]));
  end
  
  measurements = [];
  switch extension
    case '.oxyproj'
      % parse the content here
      measurements = parseoxyproj(fullfile(filepath, [filename extension]));
    case '.oxy3'
      % handling the * wildcard here
      if ismember('*', filename)
        files = dir(fullfile(filepath, [filename extension]));
        
        if isempty(files)
          error('[Artinis] There are no oxy3-files in %s. Please specify a valid filename.', filepath);
        end
        for f=files'
          measurements(end+1).file = fullfile(filepath, f.name);
          measurements(end).events = [];
          measurements(end).positions = [];
        end
      else
        % all good, single file
      end
      
    otherwise
      exts = sprintf('.%s-, ', supported_extensions{1:end-1});
      exts = sprintf('%s or .%s-', exts(1:end-3), supported_extensions{end});
      error('[Artinis] file extension not supported. please select a file in %sformat.', exts);
  end
end

if ~isempty(measurements)
  nirs_data = [];
  events = [];
  % recursively call this method for each file and save the output in a cell-array
  nf = 0;
  for f=measurements
    [p, n, e] = fileparts(f.file);
    nf = nf+1;
    if ~iscell(savename)
      savename{1} = savename;
    end
    if nargin < 2
      savename = [];
      savename{nf} = 'batch';
    elseif nargin < 3
      savename = [];
      savename{nf} = '';
    elseif numel(savename) < nf
      savename = [];
      savename{nf} = f.file;
    end
    [nirs_data{end+1} events{end+1}] = oxysoft2matlab(fullfile(p, [n e]), target, savename{nf}, verbose, f.events);
    
    % copy over 3D digitized positions if available
    for p=1:numel(f.positions)
      if ~isempty(f.positions(p).rx(1).pos)
        nirs_data{end}.positions = f.positions(p);
      end
    end
  end
  return
end

% user interface asking for saving the exported file
if nargin < 2 || strcmp(savename, 'batch') % yes, <2 - I assume that a user does not know what he does if he does not specify all input arguments
  if strcmp(target, 'homer')
    [savename,pathname] = uiputfile('.nirs', sprintf('Choose a destination for saving the export of %s to %s', filename, target));
    if savename == 0
      savename = []; % no saving desired
    else
      savename = savename(1:end-5); % will be added later
    end
  else
    [savename,pathname] = uiputfile('.mat', sprintf('Choose a destination for saving the export of %s to %s', filename, target));
    if savename == 0
      savename = []; % no saving desired
    else
      savename = savename(1:end-4); % will be added later
    end
  end
  if ~isempty(savename)
    savename = fullfile(pathname, savename);
  else
    if verbose
      fprintf('[Artinis] Save dialog aborted - export will not be saved on disk.\n');
    end
  end
elseif nargin < 3
  if verbose
    fprintf('[Artinis] command-line input detected - assuming export should not be stored on disk.\n');
  end
  savename = [];
end

%% import input file

% here, we import the relevant pieces of information. these are:
% oxy- and deoxy-values, #channels, #transs, #detector, sampling rate,
% wavelengths, optode distance, DPF (and how it was obtained)

% and, if possible, we also extract event information:
% # conditions, condition names and condition onsets
if verbose
  fprintf('[Artinis] Reading in data...');
end
switch(extension)
  case '.oxy3'
    % import to matlab function
    [rawOD,metaInfo,ADvalues] = readoxy3file(fullfile(filepath, sprintf('%s%s', filename, extension)));
    
    %     case '.txt'
    %         [rawOD, metaInfo, ADvalues] = readoxy3txtfile(fullfile(filepath, sprintf('%s%s', filename, extension)));
    %         keyboard;
    
    %case '.xml'
    %case {'.oxy', '.evt'}
  otherwise
    error('[Artinis] no support for %s-files, yet. if urgent, please send a mail to support@artinis.com.', extension);
end

%% extract data
% get all relevant datachannels by the optode template file
if verbose
  fprintf(' arranging optodes...\n')
end
[oxyOD2, metaInfo] = arrangeoxy3optodes(rawOD', metaInfo);
if verbose
  if isfield(metaInfo, 'OptodeTemplatesFile')
    fprintf('[Artinis] assuming desired optode templates in "%s...', metaInfo.OptodeTemplatesFile);
  end
end
if isnumeric(oxyOD2) && oxyOD2==-1 % if no data was found
  
  if verbose
    if isfield(metaInfo, 'OptodeTemplatesFile')
      fprintf(' but the appropriate template(s) were not found in there.\n');
    end
  end
  
  % checking for a local version (in path) of optodetemplates.xml
  metaInfo.checkLocal = true;
  [oxyOD2, metaInfo] = arrangeoxy3optodes(rawOD', metaInfo);
  if verbose
    if isfield(metaInfo, 'OptodeTemplatesFile')
      fprintf('[Artinis] assuming desired optode templates in "%s...', metaInfo.OptodeTemplatesFile);
    end
  end
  
  if isnumeric(oxyOD2) && oxyOD2==-1 % if no data was found again
    
    if verbose
      if isfield(metaInfo, 'OptodeTemplatesFile')
        fprintf(' but the appropriate template(s) were not found in there.\n');
      end
    end
    
    % let the user decide manually
    metaInfo.checkLocal = false;
    fprintf('[Artinis] no appropriate optodetemplates.xml found - please choose one manually!\n');
    [oxyOD2, metaInfo] = arrangeoxy3optodes(rawOD', metaInfo);
    if verbose
      if isfield(metaInfo, 'OptodeTemplatesFile')
        fprintf('[Artinis] assuming desired optode templates in "%s...', metaInfo.OptodeTemplatesFile);
      end
    end
    
    % okay, all previous attempts to select the correct optodetemplates.xml
    % failed and the user did not choose the correct optodetemplates.xml, so error
    if isnumeric(oxyOD2) && oxyOD2==-1
      error('[Artinis] cannot find correct optode template in optodetemplate.xml');
    end
  end
end

if verbose
  fprintf(' and I found the desired optode templates in there\n');
  fprintf('[Artinis] %d template(s) are used\n', numel(metaInfo.OptodeTemplateID));
  for i=1:numel(metaInfo.OptodeTemplateID)
    if iscell(metaInfo.OptodeTemplateID)
      fprintf('[Artinis] template#%d - ID:%d (name: %s)\n', i, metaInfo.OptodeTemplateID{i}, metaInfo.OptodeTemplateName{i});
    else
      fprintf('[Artinis] template#%d - ID:%d (name: %s)\n', i, metaInfo.OptodeTemplateID(i), metaInfo.OptodeTemplateName{i});
    end
  end
end

needsOnlyODs = any(ismember(target, {'homer', 'rawOD'}));

% concatenate all data
OD            = [];
oxyvals       = [];
dxyvals       = [];
ua            = [];
absO2Hb       = [];
absHHb        = [];
TSI           = [];
TSI_FF        = [];
nTransmitters = 0;
nReceivers    = 0;
label         = {};
RxLabel       = {};
TxLabel       = {};
TSILabel      = {};
wavelengths   = [];
distance      = [];
devices       = [];
chanPos       = [];
transPos      = [];
receiPos      = [];
transId       = [];
receiId       = [];
RxId          = [];
TxId          = [];
Rx_TxId       = [];
Rx_TxwLengths = [];
RxSubTemplateId = [];
TxSubTemplateId = [];
% note, we'll parse a lot here, in case future-me/-you needs more info
for d=1:numel(oxyOD2.Sys)
  OD            = [OD oxyOD2.Sys(d).subtemplate(:).OD];
  devices       = [devices; metaInfo.Sys(d).nRx metaInfo.Sys(d).nTx 0];
  nTransmitters = nTransmitters + metaInfo.Sys(d).nTx;
  nReceivers    = nReceivers+ metaInfo.Sys(d).nRx;
  if numel(oxyOD2.Sys)>1
    % concatenate System number here
    label         = [label strcat(sprintf('S%i-', d),[metaInfo.Sys(d).subtemplate(:).RxTx])];
    RxLabel       = [RxLabel strcat(sprintf('S%i-', d),metaInfo.Sys(d).Rx.names)];
    TxLabel       = [TxLabel strcat(sprintf('S%i-', d),metaInfo.Sys(d).Tx.names)];
  else
    label         = [label [metaInfo.Sys(d).subtemplate(:).RxTx]];
    RxLabel       = [RxLabel metaInfo.Sys(d).Rx.names];
    TxLabel       = [TxLabel metaInfo.Sys(d).Tx.names];
  end
  wavelengths   = [wavelengths metaInfo.Sys(d).Wavelength];
  distance      = [distance metaInfo.Sys(d).subtemplate(:).dist];
  chanPos       = [chanPos [metaInfo.Sys(d).subtemplate(:).pos]];
  transPos      = [transPos metaInfo.Sys(d).Tx.pos];
  receiPos      = [receiPos metaInfo.Sys(d).Rx.pos];
  
  % Below Ids are only used for giving unique numbers, thus correct for
  % multiple systems
  maxRx         = max([receiId 0]);
  maxTx         = max([transId 0]);
  transId       = [transId maxTx+[metaInfo.Sys(d).Tx.unique_id{:}]];
  receiId       = [receiId maxRx+[metaInfo.Sys(d).Rx.unique_id{:}]];
  maxRx         = max([RxId 0]);
  maxTx         = max([TxId 0]);
  RxId          = [RxId cell2mat([metaInfo.Sys(d).subtemplate(:).iRx])+maxRx];
  TxId          = [TxId cell2mat([metaInfo.Sys(d).subtemplate(:).iTx])+maxTx];
  RxSubTemplateId = [RxSubTemplateId metaInfo.Sys(d).Rx.subtemplate];
  TxSubTemplateId = [TxSubTemplateId metaInfo.Sys(d).Tx.subtemplate];
  Rx_TxId       = [Rx_TxId [metaInfo.Sys(d).subtemplate(:).cmbs]+...
    repmat([maxRx; maxTx], 1, size([metaInfo.Sys(d).subtemplate(:).cmbs], 2))];
  
  wLengths = vertcat(cell2mat([metaInfo.Sys(d).subtemplate(:).Wavelength]));
  Rx_TxwLengths = [Rx_TxwLengths; wLengths(:)];
  
  if ~needsOnlyODs || isfield(metaInfo.Sys(d).subtemplate(:), 'TSI')
    [toxyvals, tdxyvals,tTSI,tTSI_FF,tabsO2Hb,tabsHHb] = transformoxy3od([oxyOD2.Sys(d).subtemplate(:).OD],metaInfo.Sys(d),metaInfo.DPF);
    oxyvals = [oxyvals toxyvals];
    dxyvals = [dxyvals tdxyvals];
    if ~isempty(tTSI) && strcmp(target, 'oxy/dxy')
      absO2Hb = [absO2Hb tabsO2Hb'];
      absHHb  = [absHHb tabsHHb'];
      TSI     = [TSI tTSI'];
      TSI_FF  = [TSI_FF tTSI_FF'];
      for t=1:numel(metaInfo.Sys(d).subtemplate)
        if ~isempty( metaInfo.Sys(d).subtemplate(t).TSI)
          if numel(oxyOD2.Sys)>1
            % concatenate System number here
            TSILabel{end+1} = strcat(sprintf('S%i-', d), [metaInfo.Sys(d).subtemplate(t).TSI.RxTx]);
          else
            TSILabel{end+1} = metaInfo.Sys(d).subtemplate(t).TSI.RxTx;
          end
        end
      end
    end
  end
end
devices(end, 3) = metaInfo.nADC;

rawvals    = OD;
advals     = ADvalues'; % samples x channels as homer and nirs-spm
nSamples   = size(rawvals, 1);
nChannels  = numel(label); % channels == optode-pairs
nDevices   = numel(metaInfo.Device);

fs           = 1/metaInfo.SampleTime;
nWavelengths = numel(unique(wavelengths));

% differential path-length factor
DPF            = metaInfo.DPF;
DPF_correction = 'none';

if unique(length(transId))*2 ~= numel(wavelengths) && ~needsOnlyODs
  error('[Artinis] more than two laser per transmitter unit cannot be transformed to oxy- and deoxyvalues, yet. If your system had 3 wavelengths. the output of this function will be wrong! Please contact the Artinis support.');
end

%% create output structure
nirs_data = [];

switch(target)
  
  case {'rawOD', 'oxy/dxy'}
    
    if strcmp(target, 'rawOD')
      nirs_data.OD       = rawvals;
      % also change labels and channel positions respectively
      nirs_data.ODlabel = cell(1, size(rawvals, 2));
      nirs_data.ODchanPos = zeros(2, size(rawvals, 2));
      for i=1:numel(Rx_TxwLengths)/2 % TODO FIXME: hardcoded - assuming 2 wavelengths per Tx
        nirs_data.ODlabel{2*i-1} = strcat(label{i}, '@', num2str(Rx_TxwLengths(2*i-1)), 'nm');
        nirs_data.ODlabel{2*i}   = strcat(label{i}, '@', num2str(Rx_TxwLengths(2*i)), 'nm');
        
        nirs_data.ODchanPos(:, 2*i-1) = chanPos(:, i);
        nirs_data.ODchanPos(:, 2*i)   = chanPos(:, i);
      end
      
    else
      nirs_data.oxyvals  = oxyvals;
      nirs_data.dxyvals  = dxyvals;
      if ~isempty(TSI) % only if TSI data is available in the measurement
        nirs_data.absO2Hb  = absO2Hb;
        nirs_data.absHHb   = absHHb;
        nirs_data.TSI      = TSI;
        nirs_data.TSILabel = TSILabel;
        nirs_data.TSI_FF   = TSI_FF;
      end
    end
    nirs_data.time         = (0:nSamples-1)' ./ fs;
    nirs_data.wavelengths  = wLengths;
    nirs_data.DPF          = metaInfo.DPF;
    nirs_data.distance     = distance;
    nirs_data.Fs           = fs;
    nirs_data.label        = label;
    nirs_data.RxLabel      = RxLabel;
    nirs_data.TxLabel      = TxLabel;
    nirs_data.Rx_TxId      = Rx_TxId;
    nirs_data.chanPos      = chanPos;
    nirs_data.transPos     = transPos;
    nirs_data.receiPos     = receiPos;
    nirs_data.ADvalues     = advals;
    %nirs_data.metaInfo      = metaInfo;
    
    [names, onsets, durations] = parseevents();
    
    % argout
    nirs_data.events           = [];
    nirs_data.events.names     = names;
    nirs_data.events.onsets    = onsets;
    nirs_data.events.durations = durations;
    
    % argout
    events           = [];
    events.names     = names;
    events.onsets    = onsets;
    events.durations = durations;
    
    if ~isempty(savename) % save as .nirs-file
      save(savename, 'nirs_data');
      if verbose
        fprintf('[Artinis] %s data export stored as %s.\n', target, [savename '.mat']);
      end
    end
    
  case 'nirs-spm'
    
    % this is tricky, NIRS-SPM explicitly index channels starting
    % at the left in the top row, then go through the row to the
    % right and then go to the second next row, etc...
    % thus, we index according to physical position
    
    % concatenate receivers and transmitters
    optPos = [receiPos transPos];
    optIdx = [false(size(receiPos, 2), 1) receiId'; true(size(transPos, 2), 1) transId'];
    
    % sort according to y-position(pos to neg),
    % then x-position (neg to pos)
    [unused, idx] = sortrows(optPos', [-2 1]);
    optIdx = optIdx(idx, :);
    RxTx   = Rx_TxId(:, 1:2:end); % do not need the 2nd wavelength here
    
    % sort channel positions, and then everything accordingly
    [unused, idx] = sortrows(chanPos', [-2 1]);
    RxId            = RxId(:, idx);
    TxId            = TxId(:, idx);
    RxTx            = RxTx(:, idx);
    chanPos         = chanPos(:, idx);
    label           = label(:, idx);
    
    %uRxId           = unique(RxId, 'stable');
    [uA a b]        = unique(RxId,'first');
    uRxId = RxId(sort(a));
    
    %uTxId           = unique(TxId, 'stable');
    [uA a b]        = unique(TxId,'first');
    uTxId = TxId(sort(a));
    
    nirs_data.oxyData        = oxyvals(:, idx);
    nirs_data.dxyData        = dxyvals(:, idx);
    nirs_data.nch            = nChannels;
    nirs_data.fs             = fs;
    nirs_data.wavelength     = [median(Rx_TxwLengths(2:2:end)) median(Rx_TxwLengths(1:2:end))]; % is nwLengths x 1
    nirs_data.distance       = distance; % median distance, does not have to be be scalar?
    nirs_data.DPF            = repmat(DPF, size(nirs_data.wavelength)); % must be nwLengths x 1
    nirs_data.DPF_correction = DPF_correction;
    
    % extract events and names
    [names, onsets, durations] = parseevents;
    
    % argout
    events           = [];
    if ~isempty(names)
      events.names     = names;
      events.onsets    = onsets;
      events.durations = durations;
    end
    
    %% SAVING
    if ~isempty(savename) % save as .mat-file
      save([savename '_converted_data.mat'], 'nirs_data');
      if verbose
        fprintf('[Artinis] %s data export stored as %s.\n', target, [savename '_converted_data.mat']);
      end
      if ~isempty(events)
        save([savename '_multiple_conditions.mat'], 'names', 'onsets', 'durations');
        if verbose
          fprintf('[Artinis] %s event export stored as %s.\n', target, [savename '_multiple_conditions.mat']);
        end
      end
      
      % write channel configuration file
      fid = fopen([savename '_ch_config.txt'], 'w');
      % line 1: vender "Artinis Oxymon MKIII or PortaMon or PortaLite or Octamon
      fprintf(fid, 'Artinis Medical Systems Oxymon MK III\n'); % TODO FIXME
      % line 2: type, i.e. # channels, basically the subtemplate name
      templateNames = sprintf('%s & ', metaInfo.OptodeTemplateName{:});
      templateNames = templateNames(1:end-3);
      fprintf(fid, [templateNames '\n']);
      % line 3: #set(s), i.e. number of subtemplates
      if numel(metaInfo.OptodeTemplateName)> 1
        fprintf(fid, '%isets\n', size(templateNames, 1));
      else
        fprintf(fid, '1set\n');
      end
      % line 4: blank
      fprintf(fid, '\n');
      
      % from line 5 on idx of transmitter, idx of receiver (total idx in
      % data, same order as the digitizer will have!
      for i=1:nChannels
        % relabelled Rx and Tx
        Rx = find(optIdx(:, 1) == 0 & optIdx(:, 2) == RxTx(1, i));
        Tx = find(optIdx(:, 1) == 1 & optIdx(:, 2) == RxTx(2, i));
        fprintf(fid, '%i %i\n', Rx, Tx);
      end
      fclose(fid);
      
      if verbose
        fprintf('[Artinis] information about channel-configuration stored as %s.\n', [savename '_ch_config.txt']);
      end
      % a second txt file which contains the optode indices of the
      % subtemplates is not needed, the user can better do this
      % within NIRS-SPM. Alternatively, uncomment these lines:
      %             fid = fopen([savename '_ch_set.txt']);
      %
      %             for i=1:size(templateNames, 1)
      %                 fprintf('Set #%i: %i', i, find(SetRxSubTemplateId==i))%
      %             end
    end
    
    %% GRAPHICAL OUTPUT
    % as nirs-spm has no option to give visual feedback about the
    % channel positions, we do that here
    figure;
    cla;
    hold on;
    scatter(receiPos(1, :), receiPos(2, :), 180, 'r');
    scatter(transPos(1, :), transPos(2, :), 180, 'b');
    scatter(chanPos(1, :), chanPos(2, :), 180, 'k');
    done = false(size(optIdx, 1), 1);
    for i=1:nChannels
      text(chanPos(1, i), chanPos(2, i)+.5, label(i), 'HorizontalAlignment','center');
      if i<10
        text(chanPos(1, i), chanPos(2, i), num2str(i), 'HorizontalAlignment','center');
      else
        text(chanPos(1, i), chanPos(2, i), num2str(i), 'HorizontalAlignment','center');
      end
      
      r = RxTx(1, i);
      Rx = find(optIdx(:, 1) == 0 & optIdx(:, 2) == RxTx(1, i));
      if ~done(Rx)
        text(receiPos(1, r), receiPos(2, r)+.5, RxLabel(r), 'Color', 'r', 'HorizontalAlignment','center');
        if i<5
          text(receiPos(1, r), receiPos(2, r), num2str(Rx), 'Color', 'r', 'HorizontalAlignment','center');
        else
          text(receiPos(1, r), receiPos(2, r), num2str(Rx), 'Color', 'r', 'HorizontalAlignment','center');
        end
        done(Rx) = true;
      end
      
      
      
      t = RxTx(2, i);
      Tx = find(optIdx(:, 1) == 1 & optIdx(:, 2) == RxTx(2, i));
      if ~done(Tx)
        text(transPos(1, t), transPos(2, t)+.5, TxLabel(t), 'Color', 'b', 'HorizontalAlignment','center');
        if i<5
          text(transPos(1, t), transPos(2, t), num2str(Tx), 'Color', 'b', 'HorizontalAlignment','center');
        else
          text(transPos(1, t), transPos(2, t), num2str(Tx), 'Color', 'b', 'HorizontalAlignment','center');
        end
        done(Tx) = true;
      end
    end
    hold off;
    ylim(ylim+[-1 1]);
    yl = ylim;
    xlim(xlim+[-1 1]);
    xl = xlim;
    set(gca, 'XTick', xl(1):xl(end));
    set(gca, 'YTick', yl(1):yl(end));
    title('[Artinis] NIRS-SPM channel configuration and numbering');
    legend('Receiver ID', 'Transmitter ID', 'Channel ID')
    set(gcf,'color','w');
    axis on;
    grid on;
    
  case 'homer'
    
    nirs_data.t               = (0:nSamples-1)' ./ fs;
    nirs_data.d               = rawvals; % difference in OD definition ln vs log
    nirs_data.SD              = struct;
    nirs_data.SD.SrcPos       = transPos';
    nirs_data.SD.SrcPos(:, 3) = 0;
    nirs_data.SD.DetPos       = receiPos';
    nirs_data.SD.DetPos(:, 3) = 0;
    nirs_data.SD.nSrcs        = size(transPos, 2); % incorporating split fibers
    nirs_data.SD.nDets        = size(receiPos, 2); % incorporating split fibers
    
    % Homer can only handle two wavelengths, thus we need to pretend
    % each two WLs were the same across lasers
    nirs_data.SD.Lambda       = [median(Rx_TxwLengths(1:2:end)) median(Rx_TxwLengths(2:2:end))];
    
    % cannot take nChannels here, because homer wants rawdata
    nirs_data.SD.MeasList     = ones(size(rawvals, 2), 4);
    
    % re-sort so that low wavelengths preceed high wavelengths
    nirs_data.d               = 1./exp(log(10).* ...
      [rawvals(:, 1:2:end) rawvals(:, 2:2:end)]);
    
    nirs_data.SD.MeasList(1:nChannels, :) = [...
      Rx_TxId(2, 1:2:size(rawvals, 2))' ...
      Rx_TxId(1, 1:2:size(rawvals, 2))' ...
      ones(nChannels, 1) ...
      ones(nChannels, 1)];
    
    nirs_data.SD.MeasList(1+nChannels:2*nChannels, :) = [...
      Rx_TxId(2, 2:2:size(rawvals, 2))' ...
      Rx_TxId(1, 2:2:size(rawvals, 2))' ...
      ones(nChannels, 1) ...
      2*ones(nChannels, 1)];
    
    nirs_data.SD.SpatialUnit = 'cm';
    
    % this is a loop for the nonnative Matlab-speaker, doing the same
    % as just above
    %         nirs_data.SD.MeasList     = ones(size(rawvals, 2), 4);
    %         nirs_data.d = [];
    %         % HOMER wants first all lower wavelengths...
    %         j = 0;
    %         for i=1:2:size(rawvals, 2)
    %             j = j+1;
    %             nirs_data.SD.MeasList(j, 2) = Rx_TxId(1, i); % idx of det
    %             nirs_data.SD.MeasList(j, 1) = Rx_TxId(2, i); % idx of src
    %             % (i, 3) needs to be one per definition, indicating CW
    %             nirs_data.SD.MeasList(j, 4) = mod(i, 2)*-1+2; %map from [1 0] to [1 2]
    %             nirs_data.d(:, j) = 1./exp(log(10).*rawvals(:, i)); % difference in OD definition ln vs log
    %         end
    %
    %         % ... and then all higher wavelengths
    %         for i=2:2:size(rawvals, 2)
    %             j = j+1;
    %             nirs_data.SD.MeasList(j, 2) = Rx_TxId(1, i); % idx of det
    %             nirs_data.SD.MeasList(j, 1) = Rx_TxId(2, i); % idx of src
    %             % (i, 3) needs to be one per definition, indicating CW
    %             nirs_data.SD.MeasList(j, 4) = mod(i, 2)*-1+2; %map from [1 0] to [1 2]
    %             nirs_data.d(:, j) = 1./exp(log(10).*rawvals(:, i)); % difference in OD definition ln vs log
    %         end
    
    % event reading
    nirs_data.s              = zeros(nSamples, 1);
    
    nirs_data.mL             = nirs_data.SD.MeasList;
    
    
    % extract events and names
    [names, onsets, durations] = parseevents;
    
    % argout
    events           = [];
    events.names     = names;
    events.onsets    = onsets;
    events.durations = durations;
    
    for i=1:numel(unique(events.names))
      idx = ismember(events.names, events.names{i});
      nirs_data.s([events.onsets{idx}], i) = 1;
      nirs_data.CondNames{i} = events.names{i};
    end
    
    CondNames = nirs_data.CondNames;
    
    if size(advals, 2) == 0
      advals(:, 1) = 0;
    end
    nirs_data.aux = advals;
    
    %% SAVING
    if ~isempty(savename) % save as .nirs-file
      SD = nirs_data.SD;
      t = nirs_data.t;
      d = nirs_data.d;
      s = nirs_data.s;
      if size(s, 2)==0
        s(:, 1) = 0;
      end
      mL = nirs_data.mL;
      aux = nirs_data.aux;
      savename = [savename '.nirs'];
      save(savename, '-mat', 'SD','t','d', 's', 'mL', 'aux', 'CondNames');
      if verbose
        fprintf('[Artinis] %s data export stored as %s.\n', target, savename);
      end
    end
    
    % need to save this in a ".nirs" file, which is just a renamed .mat
    % save('filename.nirs', '-mat', 'SD','t','d', 's', 'ml', 'aux')
    % according to homer convention, files should be saved as
    % <subjectid>_run<n>.nirs
    % this allows homer to automatically associate individual
    % experimental sessions (here called "run") with each subject
    
    % NOTE: additonal fields could be (acc. to Homer manual Nov 29 2012, sec 6.6):
    %   tIncMan : Manual time exclusion
    %   SD.MeasListAct : Manual channel pruning
    %   CondNames : a list of condition names used for the stimulus marks
    %   s0 : The s matrix is modified as a result of stimulus mark exclusion, deletion or
    %           addition of new stimulus marks. If the s matrix is changed in any way, the original
    %           s matrix is saved in s0.
    %   userdata : Contains the information entered into the user data table in the
    %           stimGUI.
    %   procInput : See section 9.1 for details
    %   procResult : See section 9.1 for details
    % Section 9.1 is
    %   Raw data is displayed by clicking the Raw Data radio button and selecting one or both of the wavelengths in the center listbox. In the example below, raw data is displayed at 690nm.
    % (yes, that's all of Sec 9.1)
    % In contrast, Section 8.1 is:
    %   The following additional parameters will be saved in the .nirs files as a result of data processing:
    %   procInput : This variable contains the processing stream functions, their parameters and
    %           values used to process data.
    %   procResult : Contains the results of processing including optical density, concentration,
    %           optical density and concentration averages, and any other variables that are the
    %           result of executing the processing stream
    
    
  case 'nap'
    
    % NAP only supports 24 chan Oxymon
    %if iscell(metaInfo.OptodeTemplateID) || ~(metaInfo.OptodeTemplateID == 14 || metaInfo.OptodeTemplateID == 15)
    % TODO FIXME: check this, is that true?
    %error('[Artinis] the ''nap''-toolbox only supports the 24 channel oxymon system.');
    %end
    
    % re-order as [oxy1, dxy1, cmb1, oxy2, dxy2, cmb2, ...]
    nirs_data.naz_Oxy   = oxyvals;
    nirs_data.naz_Deoxy = dxyvals;
    nirs_data.naz_Total = oxyvals+dxyvals;
    
    % extract events and names
    [names, onsets, durations] = parseevents;
    
    % argout
    events           = [];
    events.names     = names;
    events.onsets    = onsets;
    events.durations = durations;
    
    Evs       = zeros(size(OD, 1), 1);
    
    if ~isempty(events) && ~isempty(events.onsets)
      Evs(cell2mat(events.onsets))= 1;
    end
    
    nirs_data.naz_Oxy(:, end+1)   = Evs;
    nirs_data.naz_Deoxy(:, end+1) = Evs;
    nirs_data.naz_Total(:, end+1) = Evs;
    
    
    %% SAVING
    if ~isempty(savename)
      naz_Oxy   = nirs_data.naz_Oxy;
      naz_Deoxy = nirs_data.naz_Deoxy;
      naz_Total = nirs_data.naz_Total;
      
      save([savename '_Oxy.mat'], 'naz_Oxy');
      if verbose
        fprintf('[Artinis] %s oxy data export stored as %s.\n', target, [savename '_Oxy.mat']);
      end
      save([savename '_Deoxy.mat'], 'naz_Deoxy');
      if verbose
        fprintf('[Artinis] %s deoxy data export stored as %s.\n', target, [savename '_Deoxy.mat']);
      end
      save([savename '_Total.mat'], 'naz_Total');
      if verbose
        fprintf('[Artinis] %s total data export stored as %s.\n', target, [savename '_Total.mat']);
      end
    end
    
    %case 'fosa'
    %case {'p3', 'potato'}
  otherwise
    error('support for %s is not yet implemented. if urgent, please send a mail to support@artinis.com.', target);
end
if verbose
  fprintf('[Artinis] please verify the correctness of your export thoroughly before using in MATLAB! Thanks!\n');
end
%% SUBFUNCTION

function [names, onsets, durations] = parseevents

if isfield(metaInfo, 'Event')
  
  %read the events from the oxy3 file
  if iscell(metaInfo.Event)
    xmlEvents = cell2mat(metaInfo.Event);
  else
    xmlEvents = metaInfo.Event(:);
  end
  
  if isfield(metaInfo, 'Description')
    [condNames,unused,evpos] = unique(metaInfo.Description);
  elseif isfield(metaInfo, 'Name')
    warning('[Artinis] extraction of event text not possible.');
    [condNames,unused,evpos] = unique(metaInfo.Name);
  else
    warning('[Artinis] extraction of event character names not possible.');
    [condNames,unused,evpos] = unique(xmlEvents);
  end
  
  
  % internally, the samples start at 0
  % this means an event at index 1 occured at sample 0
  xmlEvents = xmlEvents+1;
  
  names = condNames;
  if ~iscell(names)
    names = {names};
  end
  onsets = cell(1, numel(condNames));
  durations = cell(1, numel(condNames));
  for j=1:length(evpos)
    onsets{evpos(j)} = [onsets{evpos(j)} xmlEvents(j)];
    durations{evpos(j)} = [durations{evpos(j)} -1]; % stick function
  end
else
  names = {};
  onsets = {};
  durations = {};
end

if ~isempty(projevents)
  names = horzcat(names, projevents.names);
  onsets = horzcat(onsets, projevents.onsets);
  durations = horzcat(durations, projevents.durations);
end

end

if (verbose)
  fprintf('\n');
end
end