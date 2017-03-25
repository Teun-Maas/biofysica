function SpDe=ic_spkdenshv(Spikes,Sigma);
%
%	function [SpDe]=spkdenshv(Spikes,Sigma);
%
%	as function spkdens but normalized
% 	with 1/(sqrt(2*pi)*Sigma)
%
if Sigma==0,
   SpDe=Spikes;
else
   normfac=1/(sqrt(2*pi)*Sigma);
   SpDe=normfac*spkdens(Spikes,Sigma);
end;

