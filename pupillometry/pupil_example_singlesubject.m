
% 1) COMP_KAT_Pupil_SNR50_tr30_rb_s0n0_IFFM__02Nov2016_1501.mat ??> Adaptive procedure to look SRT
% 
% 2) COMP_KAT_SRT-16_tr20_rb_s0n0_IFFM__02Nov2016_1512.mat ??> High Intelligibility condition SRT-16dB of fix SNR.
% 
% 3) COMP_KAT_SRT+16_tr20_rb_s0n0_IFFM__02Nov2016_1526.mat ? ?> Low Intelligibility condition SRT+16dB of fix SNR.
% 
% 4) COMP_KAT_SRT_tr20_rb_s0n0_IFFM__02Nov2016_1519.mat ? ?> Medium Intelligibility condition SRTdB of fix SNR. 
close all
clearvars;

cd('/Users/marcw/DATA/Sebastian Ausili/pupil');

% fnames = {'COMP_KAT_Pupil_SNR50_tr30_rb_s0n0_IFFM__02Nov2016_1501.mat',...
% 'COMP_KAT_SRT-16_tr20_rb_s0n0_IFFM__02Nov2016_1512.mat',...
% 'COMP_KAT_SRT+16_tr20_rb_s0n0_IFFM__02Nov2016_1526.mat',...
% 'COMP_KAT_SRT_tr20_rb_s0n0_IFFM__02Nov2016_1519.mat'};
% 
% 
% cond	= {'adapt','high','low','medium'};
% SNR		= [0 +16 -16 0];

fnames = {'COMP_Martijn_Pupil_tr30_rb_s0n0_IFFM__07Nov2016_1622.mat',...
'COMP_Martijn_Pupil_tr20_rb_s0n0_IFFM__07Nov2016_1649.mat',...
'COMP_Martijn_Pupil_tr20_rb_s0n0_IFFM__07Nov2016_1637.mat',...
'COMP_Martijn_Pupil_tr20_rb_s0n0_IFFM__07Nov2016_1643.mat',...
'COMP_Martijn_Pupil_tr20_rb_s0n0_IFFM__07Nov2016_1631.mat'};
 
 
cond    = {'adapt','silence','high','medium','low'};
SNR     = [0 +100 +20 -20 0];

	 [data_out,rec] = pupil_preprocess(fnames{4},'disp',true);
return
	%%
nfiles = numel(fnames);
for ii = 1:nfiles
	fname = fnames{ii};
	data_out(ii) = pupil_preprocess(fname,'disp',false);
% 	pause
end


%%
col = lines(nfiles);

MPD = [data_out.MPD];
PPD = [data_out.PPD];
L	= [data_out.latency];

close all
subplot(222)
plot(SNR,MPD,'ko','MarkerFaceColor','w')
axis square;
box off
xlim([-20 20]);
% ylim([-0.05 0.35]);
set(gca,'TickDir','out',...
	'XTick',-16:4:16,'YTick',0:0.1:0.3);
xlabel('SNR re SRT (dB)');
ylabel('mean pupil dilation (mm)');

subplot(224)
plot(SNR,L,'ko','MarkerFaceColor','w')
axis square;
box off
xlim([-20 20]);
% ylim([1 3]);
set(gca,'TickDir','out',...
	'XTick',-16:4:16,'YTick',1:0.5:3);
xlabel('SNR re SRT (dB)');
ylabel('latency peak pupil dilation (s)');

subplot(223)
plot(SNR,PPD,'ko','MarkerFaceColor','w')
axis square;
box off
xlim([-20 20]);
% ylim([-0.05 0.45]);
set(gca,'TickDir','out',...
	'XTick',-16:4:16,'YTick',0:0.1:0.4);
xlabel('SNR re SRT (dB)');
ylabel('peak pupil dilation (mm)');


subplot(221)
cla
	hold on
for ii = 1:nfiles
	time = data_out(ii).time;
	avg_data = data_out(ii).avg;
	se_data = data_out(ii).se;
	
% 	plot(time,avg_data);
	h(ii) = errorpatch(time,avg_data,se_data,col(ii,:));

end
xlim ([-1.5 5.2]);
	horline(0,'k-');
	verline([-1 0 5],'k:');
	axis square;
	xlabel('Time (s)');
	ylabel('Pupil diameter (mm)');
	box off;
	set(gca,...
		'YTick',0:0.1:0.3,...
		'XTick',-1:1:5,...
		'TickDir','out');
	
legend(h,cond,'Location','SW')

savegraph(mfilename,'eps');