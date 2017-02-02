function omega = bf_mode(x)
% OMEGA = BF_MODE(X)
%
% Determine mode, or most frequent value, OMEGA for vector X, which are
% samples from a continuous distribution.
%
% From Matlab's MODE:
% "Applying the mode function to a sample from that distribution is
% unlikely to provide a good estimate of the peak; it would be better to
% compute a histogram or density estimate and calculate the peak of that
% estimate."
%
% see also KSDENSITY, MODE

[f,xi]		= ksdensity(x);
[~,indx]	= max(f);
omega		= xi(indx);