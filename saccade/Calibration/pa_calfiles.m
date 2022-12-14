function pa_calfiles(DatFiles,CalFile)
% Calibrate all raw data contained in one or several files
%
% CALIBRATE(<DATFILES>,<CALFILE>)
%
%  Function to calibrate raw data into azimuth/elevation angles. These
%  angles are stored in hv-files.
%
%
%
%       DATFILES:       Raw data. 1 or more files in rows.
%                   eg. ['XX-XX-2000-01-0101<.dat>'; ...
%                        'XX-XX-2000-01-0102<.dat>'; ...
%                        'XX-XX-2000-01-0103<.dat>'];
%                       default: all dat-files in current directory
%
%       CALFILE:        Neural network file.
%                   eg. 'XX-XX-2000-01-01<.net>';
%                       default: user input
%
%
%       output:         hv-files named according to the DatFiles.
%                   eg. 'XX-XX-2000-01-0101.hv' ... etc.
%
%  See also ULTRADET, TRAINCAL
%
%  Author: Marcus
%  Date: 11-04-07


%% Initialization
if nargin<1
    d                        = dir('*.dat'); % Changes dir('*.*at') to dir('*.dat') on 18 aug 2022 MW
    DatFiles                = char(d.name);
end
if nargin<2
    d                        = dir('*.net');
    CalFile                = char(d.name);
    CalFile                 = pa_fcheckexist(CalFile,'*.net');
end

%
warning('MATLAB:pa_calfiles:channeldimension','Is your calibration correct? Check pa_calfiles');
% Original dimensions (for sphere??)
% Hchan                       = 3;
% Vchan                       = 1;
% Fchan                       = 2;

Hchan                       = 1;
Vchan                       = 2;
Fchan                       = 3;
S                           = load(CalFile,'-mat');



%% Calibrate all DATfiles
for i                       = 2:size(DatFiles,1)
    % Loading file
    fname                   = pa_fcheckext(DatFiles(i,:),'.dat');
%     fname                   = pa_fcheckexist(DatFiles(i,:),'*.*at');
    csvname                 = pa_fcheckext(DatFiles(i,:),'.csv');
    [expinfo,chaninfo]      = pa_readcsv(csvname);
    nchan                   = expinfo(1,8);
    nsample                 = chaninfo(1,6);
    DAT                     = pa_loaddat(fname,nchan,nsample);
    H                       = squeeze(DAT(:,:,Hchan));
    H                       = H(:);
    V                       = squeeze(DAT(:,:,Vchan));
    V                       = V(:);
    F                       = squeeze(DAT(:,:,Fchan));
    F                       = F(:);
    DAT                     = [H V F]';
	
    [AZ,EL]                 = pa_calib(DAT,S);

    % Saving calibrated data
    fname                   = pa_fcheckext(DatFiles(i,:),'.hv');
    fid                     = fopen(fname,'w','l');
    AZEL                    = [AZ;EL];
    fwrite(fid,AZEL,'float');
    fclose(fid);
end


