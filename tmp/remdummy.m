function Stim = remdummy(Stim)
% STIM = REMDUMMY(STIM)
%
% Remove dummy movements from STIM matrix
%

% (c) Marc van Wanrooij 2008
sel     = Stim(:,3)==0;
Stim    = Stim(~sel,:);
t       = unique(Stim(:,1));
for i   = 1:length(t)
    sel = Stim(:,1) == t(i);
    Stim(sel,1) = i;
end