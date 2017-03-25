% stim=fgetstim(fpLOG,TrialNr,Ntarget)
%   read stimuli from log file. TrialNr is a vector of 
%   trial numbers to be read. Ntarget is the number of
%   targets in a trial. The returned matrix stim is 
%   sorted on trial number. See Index for layout of 
%   stim matrix.
%
%   Jeroen Goossens

function stim=ic_fgetstim(fpLOG,TrialNr,Ntarget)

TrialNr = sort(TrialNr);

frewind(fpLOG);
i=1;
while ~feof(fpLOG), 
  LineNr = sprintf('#%d',TrialNr(i));
  ident=fscanf(fpLOG,'%s',1);
  if strcmp(ident,LineNr),
    % read line of targets 
    t1=[];
    t2=[TrialNr(i),Ntarget];
    for j=1:Ntarget,
      t1=fscanf(fpLOG,'%f',7)';
      t2=[t2, t1(1:3), pa_rphi2azel([t1(2) t1(3)]), t1(4:7)]; %PH PATCH
    end;
    stim(i,:)=t2;
    if i==length(TrialNr), 
      return; 
     end;
    i=i+1;
  end;
end;