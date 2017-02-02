function spheremat2hoopdat(fname)
% SPHEREMAT2HOOPDAT(FNAME)
%
% Convert the SPHERE MAT-file structure to the DAT-file structure of the
% HOOP
%
% See also SPHEREMAT2HOOPCSV

if nargin<1
	fname = [];
end

%% Initialization
ext             = '.mat';
matfile         = fcheckexist(fname,ext);
matfile         = fcheckext(matfile,ext);


ext             = '.dat';
datfile         = fcheckext(matfile,ext);

%% Load data
load(matfile)
ntrials = numel(data);

%% Create dat-file
fid = fopen(datfile,'w','l');
for ii = 1:ntrials
	d		= data(ii).raw;
	d = d(:,[2 3 1 4:8]); % the configuration of the SPHERE coils is different from HOOP coils
	fwrite(fid,d,'float');
end
fclose(fid);


