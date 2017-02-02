function cfg1 = addfields(cfg1,cfg2)
% STRUCT1 = ADDFIELDS(STRUCT1,STRUCT2)
%
% Add fields of one structure STRUCT2 to another structure STRUCT1
str		= fieldnames(cfg2);
for ii	= 1:numel(str)
	cfg1.(str{ii}) = cfg2.(str{ii});
end