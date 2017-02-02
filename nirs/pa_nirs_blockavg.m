function [MU,baseline] = pa_nirs_blockavg(nirs,chanSig,mod)
fs			= nirs.Fs;
fd			= nirs.fsdown;
onSample	= ceil([nirs.event.sample]*fd/fs); % onset and offset of stimulus
offSample	= onSample(2:2:end); % offset
onSample	= onSample(1:2:end); % onset
% stim		= {nirs.event.stim}; % stimulus modality - CHECK HARDWARE WHETHER THIS IS CORRECT
stim		= nirs.event.stim; % stimulus modality - CHECK HARDWARE WHETHER THIS IS CORRECT
if numel(stim)<3
	stim = {nirs.event.stim};
end

selOn		= strcmp(stim,mod);
selOff		= selOn(2:2:end);
selOn		= selOn(1:2:end);
Aon			= onSample(selOn);
Aoff		= offSample(selOff);

mx			= min((Aoff - Aon)+1)+200;
nStim		= numel(Aon);
MU = NaN(nStim,mx);
for stmIdx = 1:nStim
	idx				= Aon(stmIdx)-100:Aoff(stmIdx)+100; % extra 100 samples before and after
	idx				= idx(1:mx);
	MU(stmIdx,:)	= chanSig(idx);
end
baseline = nanmean(MU(:,1:100),2);
MU = bsxfun(@minus,MU,baseline); % remove the 100th sample, to set y-origin to stimulus onset

