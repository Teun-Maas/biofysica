% function pa_johnic_predict
home;
clear all;
close all;
clc

warning off;
%% IC
% d = 'E:\DATA\John_IC\predictions2';
% cd(d);


%% STRF
strfFiles = 'E:\MATLAB\PANDA\Spike\IC\strfbr0812.mat';
if ~exist(strfFiles,'file')
	strfFiles = 'C:\MATLAB\PANDA\Spike\IC\strfbr0812.mat';
end
if ~exist(strfFiles,'file')
	strfFiles = 'D:\MATLAB\PANDA\Spike\IC\strfbr0812.mat';
end
load(strfFiles);

%% Vocalizations
sndFiles = 'E:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
if ~exist(sndFiles,'file')
	sndFiles = 'C:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
end
if ~exist(sndFiles,'file')
	sndFiles = 'D:\MATLAB\PANDA\Spike\STRF\timetraces2.mat';
end
load(sndFiles)

%%
fname = 'E:\MATLAB\PANDA\Spike\IC\BR08\BR0814.AAP';
if ~exist(fname,'file')
fname = 'C:\MATLAB\PANDA\Spike\IC\BR08\BR0814.AAP';
end
if ~exist(fname,'file')
	sndFiles = 'D:\MATLAB\PANDA\Spike\IC\BR08\BR0814.AAP';
end

ic_predictcall(strf, fname, h, 1);