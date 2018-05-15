% INIT_BIOFYSICA - run commands for initialization of the biophysics toolbox
% call this from your startup.m by adding two lines
%    addpath('/Users/gunter/Documents/Matlab/biofysica');
%    init_biofysica;

% init files are stored in biofysica/init.d and run in alphabetical
% order. Architecture specific files are in
% biofysica/init.d/{glnx64,maci64,win64}

[BIOFYSICA_ROOT,~,~]=fileparts(mfilename('fullpath'));
INIT_ARCH=computer('arch');
INIT_DIR=strcat(BIOFYSICA_ROOT,'/init.d');
INIT_ARCH_DIR=strcat(BIOFYSICA_ROOT,'/init.d/',INIT_ARCH);

INIT_CWD=cd;

% Run generic commands
cd(INIT_DIR);
INIT_FILES=dir('*.m');
for n=1:size(INIT_FILES)
    [~,INIT_CMD,~]=fileparts(INIT_FILES(n).name);
    eval(INIT_CMD);
end

% Run architecture specficic commands
cd(INIT_ARCH_DIR);
INIT_FILES=dir('*.m');
for n=1:size(INIT_FILES)
    [~,INIT_CMD,~]=fileparts(INIT_FILES(n).name);
    eval(INIT_CMD);
end

cd(INIT_CWD);
clear n INIT_*
