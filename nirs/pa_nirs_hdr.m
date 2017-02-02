function [hemo, X] = pa_nirs_hdr(nirs)
% Generates predicted hemodynamic response
fs		= nirs.fsample;
fd		= nirs.fsdown;
R		= nirs.processed(2,:);
on		= ceil([nirs.event.sample]*fd/fs);
off		= on(2:2:end);
on		= on(1:2:end);
N		= length(R);
X       = zeros(1, N);
for ii	= 1:length(on)
	X(on(ii):off(ii)) = 1;
end

plot(X,'k-')
% beta0 = Hb/10; % gain, offset
% beta  = nlinfit(X,R,@hemofunction,beta0);

hemo = pa_nirs_hdrfunction(1,X);
hold on
plot(hemo,'r-');