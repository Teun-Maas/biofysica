% INIT_BIOFYSICA - run commands for initialization of the biophysics toolbox
% call this from your startup.m by adding two lines
%    addpath('/Users/gunter/Documents/Matlab/biofysica');
%    init_biofysica;

% init files are stored in biofysica/init.d and run in alphabetical
% order. Architecture specific files are in
% biofysica/init.d/{glnx64,maci64,win64}

[BIOFYSICA_ROOT,~,~]=fileparts(mfilename('fullpath'));
ARCH=computer('arch');
INIT_DIR=strcat(BIOFYSICA_ROOT,'/init.d');
INIT_ARCH_DIR=strcat(BIOFYSICA_ROOT,'/init.d/',ARCH);

INIT_CWD=cd;

% Run generic commands
cd(INIT_DIR);
FILES=dir('*.m');
for n=1:size(FILES)
    [~,CMD,~]=fileparts(FILES(n).name);
    eval(CMD);
end

% Run architecture specficic commands
cd(INIT_ARCH_DIR);
FILES=dir('*.m');
for n=1:size(FILES)
    [~,CMD,~]=fileparts(FILES(n).name);
    eval(CMD);
end

cd(INIT_CWD);
clear INIT_CWD INIT_DIR INIT_ARCH_DIR FILES CMD ARCH

