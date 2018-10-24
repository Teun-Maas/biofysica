function result=lsl_pupil_pyrepr2mat(srepr)
    % LSL_PUPIL_PYREPR2MAT - convert a pupil labs PYTHON REPR string to a matlab structure.
    % A REPR string can be captured using the pupil_lsl plugin in PupilCapture and
    % selecting one of the Pupil Python Representation streams.
    % The python srepr is recursively converted to a matlab structure.
    % If SREPR is an array of PYTHON REPR strings, the result is an array
    % of matlab structure. You can convert this to a structure of arrays
    % using the AOS2SOA function.

    % GW/20180606

    global biofpy

    if ischar(srepr)
        dictvar=biofpy.eval_srepr(srepr);
        result=py_cast_recursive(dictvar);
    elseif iscell(srepr)
        result=convert_multi(srepr);
    else
        error('oops, cannot handle data of class %s',class(srepr));
    end
end

function r=convert_multi(s)
    [m,n]=size(s);
    r=cell(m,n);
    for ii=1:m
       for jj=1:n
          r{ii,jj}=lsl_pupil_pyrepr2mat(s{ii,jj});
       end
    end
    r=[r{:}];  % convert the cell array in r to a structure array
end

function result=py_cast_recursive(var)
    
    vartype=class(var);
    switch vartype
        
        case {'py.int','py.long', 'py.array.array'}
            result=double(var);
            
        case 'py.bytes'
            result=int8(var);
            
        case {'py.str','py.unicode'}
            result=string(var);
            
        case 'py.dict'
            result=cast_dict(var);

        case 'py.tuple'
            result=cast_tuple(var);
            
        case 'py.list'
            result=cast_list(var);
            
        otherwise
            if strcmp(vartype(1:3),'py.')
                warning('I don''t know this python class (yet): %s', vartype);
            end
            result=var;      
    end
end

function result=cast_dict(pyvar)
    var=struct(pyvar);
    fieldnames=fields(var);
    S.type='.';
    result=struct;
    n=numel(fieldnames);
    for ii=1:n
       S.subs=fieldnames{ii};
       pyvalue=var.(fieldnames{ii});
       result.(fieldnames{ii})=py_cast_recursive(pyvalue);
    end
end

function result=cast_tuple(pyvar)
    var=cell(pyvar);
    n=numel(var);
    result=cell(1,n);
    for ii=1:n
        pyvalue=var{ii};
        matvalue=py_cast_recursive(pyvalue);
        result{ii}=matvalue;
    end
    
    try
       result=cell2mat(result);
    catch
       % nice try, conversion not possible, proceed.
    end
end
  
function result=cast_list(pyvar)
    result=cast_tuple(pyvar); 
end





