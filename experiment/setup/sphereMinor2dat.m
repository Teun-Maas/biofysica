function sphereMinor2dat(fname)
% SPHERE2HOOPDAT(FNAME)
%
% Convert the SPHERE-file structure to the DAT-file structure of the
% HOOP
%
% See also SPHEREMAT2HOOPCSV

if nargin<1
	fname = [];
end

%% Initialization
ext             = '.sphere';
matfile         = fcheckexist(fname,ext);
matfile         = fcheckext(matfile,ext);


ext             = '.dat';
datfile         = fcheckext(matfile,ext);

%% Load data
load(matfile,'-mat')
ntrials = numel(data);

%% Create dat-file
fid = fopen(datfile,'w','l');
for ii = 1:ntrials
	d		= data(ii).raw*1700;
% 	d = d(:,[2 3 1 4:8]); % the configuration of the SPHERE coils is different from HOOP coils
	fwrite(fid,d,'float');
end
fclose(fid);


