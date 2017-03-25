
%
% function [Stim, Spk] = apestim(DatFile)
%
% Get stim and spike-event matrix.
%
% .. Paul..


% Huib Versnel/John van Opstal/Marcel Zwiers
% Copied 2012 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com


function [Stim, Spk] = ic_apestim(DatFile)

LogFile = pa_fcheckext(DatFile,'.log');

disp (['-> Loading ' DatFile]);

%%%%%%%%%%%%%%%%%%%%%% load log data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LogFile = [DatFile '.log'];
FpLog = fopen(LogFile,'rt','l');
if FpLog==-1,
  disp(['   Error opening file ' LogFile]);
  return;
end;
% get header information
Nchannel = ic_fgetpar(FpLog,'Sig.Nchannel','%d');
Nsample  = ic_fgetpar(FpLog,'Sig.Nsample','%d');
Fsample  = ic_fgetpar(FpLog,'Sig.Fsample','%d');
Tsample  = 1000/Fsample;
Ntest    = ic_fgetpar(FpLog,'Exp.Ntest','%d');
Npresent = ic_fgetpar(FpLog,'Exp.Npresent','%d');
Ntarget  = ic_fgetpar(FpLog,'Exp.Ntarget','%d');
Ntrial   = Ntest*Npresent;

% read trials
Stim  = ic_fgetstim(FpLog,[1:Ntrial],Ntarget);
Nstim = size(Stim,1);
fclose(FpLog);

%%%%%%%%%%%%%%%%%%%%%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DatFile = pa_fcheckext(DatFile,'.AAP');
% if Nchannel==6,
%   [I,S1,S2]=fgetsig(DatFile,Nchannel,Nsample,[1,2,5],[1:Ntrial],1,Nsample);
% end;
% if Nchannel==4,
  [S1,S2]=ic_fgetsig(DatFile,Nchannel,Nsample,[2,4],[1:Ntrial],1,Nsample);
% end;
Spk = ic_spkbins(S1,S2);




