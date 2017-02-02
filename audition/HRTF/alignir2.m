%
% function ir2 = alignir2 (ir)
%
% Align the impulseresponses by using cross correlations.
% The average impulse response is taken as a reference.
% The rest is aligned such that the maximum of the 
% correlation falls at zero.
%
% Use this after align, so that at least a rough alignment
% is already done.
%  
% .. Paul ...

function IR2 = pa_alignir2 (IR)

IR0 = IR;
IR  = IR0(30:100,:);

[NFFT,Nref] = size(IR);

% find the shift dn with respect to a neighboring 
% impulse response 

DN = zeros(Nref,1);
for i=2:Nref
  irref = IR(:,i-1);
  ir    = IR(:,i);

  [~,n] = max(xcorr(ir,irref, 'unbiased'));
  dn = n - NFFT;
  DN(i) = DN(i-1) + dn;
end;

% align all with the cumulative shift DN

IR2 = IR0;
for i=2:Nref
  IR2(:,i) = pa_shift(IR0(:,i),DN(i));
end;

