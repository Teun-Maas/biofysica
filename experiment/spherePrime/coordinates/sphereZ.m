function trial = sphereZ(trial,cfg)
% TRIAL = SPHEREZ(TRIAL,CFG)
% Convert azimuth, elevation into Z values for sphere
%
% Note that X,Y values do not exist for the SPHERE, only the HOOP
%
% Should be made more accesible
% e.g. not based on structures
%

for trlIdx	= 1:cfg.ntrials
	for stmIdx	= 1:trial(trlIdx).nstim
		X			= trial(trlIdx).stim(stmIdx).X;
		if ~isempty(X)
			ZI										= cfg.interpolant(trial(trlIdx).stim(stmIdx).azimuth,trial(trlIdx).stim(stmIdx).elevation);
			trial(trlIdx).stim(stmIdx).Z			= ZI;
			trial(trlIdx).stim(stmIdx).azimuth		= cfg.lookup(ZI+1,5);
			trial(trlIdx).stim(stmIdx).elevation	= cfg.lookup(ZI+1,6);
			
		end
	end
end
