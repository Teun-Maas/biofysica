function runExperimentTTL(handles)

% get trial information
handles = setupShow(handles);
handles		= gethandles(handles);
handles		= gettrial(handles);
createdir(handles);

%rc = run_pupil; % run 'pupil_remote_control'

tstart		= tic;

%%
trialClean(struct([]),handles.cfg);
for trlIdx	= 1:handles.cfg.ntrials
	t				= tic;

	
	%% Random intertrial interval
	r				= rndval(handles.cfg.ITI(1),handles.cfg.ITI(2),1);
	handles.trial(trlIdx).ITI		= r; % ms
	pause(r/1000);
	
	%% Trial number
	handles.cfg.trial		= trlIdx;
	disp(['Trial: ' num2str(trlIdx)])
	
	%% Trial
	stim			= handles.trial(trlIdx).stim;
	trialClean(stim,handles.cfg);
	
	[stim,handles.cfg,stimttl]	= trialSetupTTL(handles.cfg,stim);
	trialRunTTL(handles.cfg,stim,stimttl);
	handles.data		= trialSave(handles.cfg,handles.trial,handles.data);
	handles				= trialShow(handles);
	trialClean(stim,handles.cfg);
	
	%% End
    
     % rc.stop_recording;
    
	dur		= toc(t); %#ok<NASGU>
	% save timing
	[~,fname,~]		= fileparts(handles.cfg.fname); % remove extension
	fname			= [fname '-' num2str(handles.cfg.trial,'%04u')]; %#ok<AGROW>
	fname			= fcheckext(fname,'sphere');
	fname			= fullfile([handles.cfg.dname filesep 'trial' filesep],fname);
	save(fname,'dur','-append');
    
    
	toc(t)
end
handles.cfg.duration	= duration(0,0,toc(tstart)); % add duration of experiment

%% Finish
str = {	'DATA ARE SAVED PER TRIAL. SEE ALSO SPHERETRIAL2COMPLETE',...
	'______________________________________________________',...
	' ',...
	['Experiment was completed in: ' char(handles.cfg.duration) ' hours:minutes:seconds'],...
	['Files were saved to ' handles.cfg.dname filesep 'trial'],...
	};
str = char(str);
disp(str);
cd(handles.cfg.dname)

endBlockTTL(handles.cfg);