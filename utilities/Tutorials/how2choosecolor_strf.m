close all
clear all
clc

pa_genripple(4,0,100,1000,500,'display',1);

%%
cd('E:\DATA\Test'); % go to data-directory
% load('test.mat')
% whos
% spike = spikep;
% spike = rmfield(spike,'spikewave');
% save examplecell spike;
load('examplecell'); % load data

data = pa_spk_ripple2strf(spike,'shift',1);




pa_spk_plotstrf(data);
figure(1)
col = pa_statcolor(64,[],[],[],'def',6,'disp',false);
colormap(col);

figure(200)
col = pa_statcolor(64,[],[],[],'def',8,'disp',false);
set(gca,'XTick',0:20:100,'YTick',0:0.5:2.5);
	
print('-depsc','-painter',[mfilename 'rgb']);

colormap(col);
print('-depsc','-painter',[mfilename 'hcl']);

