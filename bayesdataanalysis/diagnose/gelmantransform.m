function x = gelmantransform(x)
% X = GELMANTRANSFORM(X)
%
% X is a MCMC structure
%
% Gelman and Rubin diagnostic assumes a normal distribution. To
% improve the normal approximation, variables on [0, Inf) are log
% transformed, and variables on [0,1] are logit-transformed.

parNames		= fieldnames(x);

if nvar(x) == 1 % for one variable/parameters
	z = x.(parNames);
	if min(z) > 0
		if(max(z) < 1) % on [0,1]
			y = log(z/(1-z)); % logit transform
		else % on [0, Inf]
			y = log(z); % log transform
		end
		x.(parNames) = y;
	end
	% do not transform if any<0
else % for multiple variables
	for i = 1:nvar(x)
		z = x.(parNames{i});
		if min(z) > 0
			if max(z) < 1
				y = log(z/(1 - z));
			else
				y = log(z);
			end
			x.(parNames{i}) = y;
		end
	end
end