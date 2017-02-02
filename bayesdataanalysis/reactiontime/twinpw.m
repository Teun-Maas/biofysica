function PW = twinpw(T,l1,l2,ga,varargin)
% PW = TWINPW(T,L1,L2,GA)
%
% Probability of warning PW of the time-window-of-integration (TWIN)
% model, given stimulus onset asynchrony T, rate of the target L1, and rate
% of the distractor (in focused attention tasks) or target (in redundant
% target tasks) L2, and GA.
%
% For model description see e.g. ???
%
% See also TWINPI

task = keyval('task',varargin,'FA'); % task = focused attention (FA) or redundant target (RT)
n	= numel(T);
PW	= NaN(n,1);


switch task
	case 'FA'
		for tIdx = 1:n
			t = T(tIdx);
			if (t+ga)<0
				PW(tIdx) = 1-(l1./(l1+l2).*exp(l2*(t+ga)));
			elseif (t+ga)>=0
				PW(tIdx) = (l2./(l1+l2).*exp(-l1*(t+ga)));
			end
		end
	case 'RT'
		disp('Not yet implemented');
end
