
function ir2 = pa_alignir(ir, ndelay)

% constants

[NFFT,Nir]   = size(ir);

% find the peaks

[~,nactual] = max (ir);
dN = ndelay - nactual;


IsNeg = find(dN<0);
if ~isempty(IsNeg)
  dN(IsNeg) = dN(IsNeg) + NFFT;
end;

ir2 = ir;
for i=1:Nir
  ir2(:,i) = pa_shift (ir(:,i), dN(i));
end;
