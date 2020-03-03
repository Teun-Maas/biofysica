function y = replacexy(x,new,old)
% Y = REPLACEXY(X,NEW,OLD)
%
% Replace OLD values in X with NEW values for numerical arrays.
%
% For structures will replace OLD fieldnames with NEW fieldnames.

if isnumeric(x)
	n		= numel(new);
	y		= x;
	for ii	= 1:n
		y(x == old(ii)) = new(ii);
	end
	
elseif isstruct(x)
	%% For structures
	y				= x;
	if iscell(new)
		% 	disp('rename fields in structure')
		nFields			= numel(new);
		for ii = 1:nFields
			oldField		= old{ii};
			newField		= new{ii};
			[y.(newField)]	= y.(oldField);
			y				= rmfield(y,oldField);
		end
	elseif ischar(new)
		% 		disp('rename only 1 field in structure')
		oldField		= old;
		newField		= new;
		[y.(newField)]	= y.(oldField);
		y				= rmfield(y,oldField);
	end
end
