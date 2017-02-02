% Example for BINORATE

%% Data from Doing Bayesian Data Analysis
DBDA = '/Users/marcw/Gitlab/thirdparty/DBDA2Eprograms';
addpath(genpath(DBDA));



%% Read the data
fname = which('BattingAverage.csv');
myData	= csvread(fname,1,2);

% %-------------------------------------------------------------------------------
% % Generate the MCMC chain:
% startTime = proc.time()
% mcmcCoda = genMCMC( data=myData ,
%                     zName="Hits", NName="AtBats", sName="Player", cName="PriPos",
%                     numSavedSteps=11000 , saveName=fileNameRoot ,
%                     thinSteps=20 )
% stopTime = proc.time()
% elapsedTime = stopTime - startTime
% show(elapsedTime)
% %-------------------------------------------------------------------------------

