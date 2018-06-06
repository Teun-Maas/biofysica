function result=lsl_pupil_srepr2mat(srepr)
    % LSL_PUPIL_SREPR2MAT - convert a pupil labs SREPR string to a matlab structure.
    % An SREPR string can be captured using the pupil_lsl plugin in PupilCapture and
    % selecting one of the Pupil Python Representation streams.
    % The python srepr is recursively converted to a matlab structure.

    % GW/20180606

    global biofpy
    dictvar=biofpy.eval_srepr(srepr);
    result=py_cast_recursive(dictvar);
 
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





