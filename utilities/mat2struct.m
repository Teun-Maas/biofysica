function r = mat2struct(m, fieldlist)
% mat2struct Create a structure from a matrix
% mat2struct(M,FIELDS) returns a structure with fieldnames
% from the cell array FIELDS filled with the rows from matrix M
% if FIELDS is a struct, it will be used as a template for the result.
    assert(ismatrix(m));    

    if isstruct(fieldlist)
        % use it as a template
        fieldlist = fieldnames(fieldlist); 
    end
    assert(iscell(fieldlist));
    
    mrow = size(m,1);
    nfield = length(fieldlist);
    assert(mrow>=nfield);
    
    if mrow > nfield
        warning('discarding %d extra columns in input matrix',mrow-nfield);
    end
    r = struct();
    for ii = 1:nfield
        r.(fieldlist{ii}) = m(ii,:);
    end
end