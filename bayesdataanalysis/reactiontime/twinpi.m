function PI = twinpi(T,l1,l2,o,varargin)
% PI = TWINPI(T,L1,L2,O)
%
% Probability of integration PI of the time-window-of-integration (TWIN)
% model, given stimulus onset asynchrony T, rate of the target L1, and rate
% of the distractor (in focused attention tasks) or target (in redundant
% target tasks) L2, and width of the time window O.
%
% For model description see e.g. ???
%
% See also TWINPW

task = keyval('task',varargin,'FA'); % task = focused attention (FA) or redundant target (RT)
n	= numel(T);
PI	= NaN(n,1);

switch task
	case 'FA'
		for tIdx = 1:n
			t = T(tIdx);
			if t<(t+o) && (t+o)<0
				PI(tIdx) = l1./(l1+l2).*exp(l2.*t).*(-1+exp(l2.*o));
			elseif t<=0 && 0<=(t+o)
				PI(tIdx) = 1./(l1+l2).*(l2*(1-exp(-l1.*(o+t)))+l1*(1-exp(l2.*t)));
			elseif 0<t && t<(t+o)
				PI(tIdx) = l2./(l1+l2).*(exp(-l1*t)-exp(-l1*(o+t)));
			end
		end
	case 'RT'
		disp('Not yet implemented');
end
