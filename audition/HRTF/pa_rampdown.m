function Sig = pa_rampdown(Sig, nn)
%
% function SIG = RAMPDOWN (SIG, [n1 n2])
%
% Modulate everything between sample/row n1 and n2
% with a squared sine (going down from 1 to 0).
% Everything after n2 is set to zero.


% .. Drs. P ...
% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

[M,N] = size(Sig);

n1 = nn(1);
n2 = nn(2);

dn = n2-n1+1;

if M<dn || ~pa_inrange(n1,[1 M]) || ~pa_inrange(n2,[1 M]) || n2<=n1
	disp ('Error (RAMPUP): illegal parameters');
else
	x		= (0:(dn-1))/(dn-1);
	Env	= ( sin(0.5*pi*x) ).^2;
	Env	= Env(end:-1:1);
	
	for i=1:N
		Sig(n1:n2,i) = Env' .* Sig(n1:n2,i);
	end;
	if (n2 < M)
		Sig(n2:end,:) = 0;
	end;
	
end;

