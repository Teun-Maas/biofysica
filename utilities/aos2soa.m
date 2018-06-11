function soa=aos2soa(aos)
    % B=AOS2SOA(A) - recursive conversion of an Array Of Structures in A to
    % a Structure Of Arrays in B.
    
    % GW/20180611
    soa=recursive_copy_struct(aos);
end


function dst_data=recursive_copy_struct(s,srcpath)
    if ~exist('srcpath','var') || isempty(srcpath)
        srcpath={};
    end
    % use s(1).srcpath.... as a template for traversing the data
    src_data=getfield(s,{1},srcpath{:});
    dst_data=struct();
    
    assert(isstruct(src_data(1)),'expect a struct in arg1')
        
    subfields=fieldnames(src_data(1));
    n=numel(subfields);
    for i=1:n
        field=subfields{i};
        current_path=srcpath;
        current_path{end+1}=field;
        if isstruct(src_data(1).(field))
            tmp=recursive_copy_struct(s,current_path);
            dst_data=setfield(dst_data,current_path{:},tmp);
        else
            tmp=arrayfun(@(x) getfield(x,current_path{:}),s,'UniformOutput',false);
            % tmp is a cell array now
            try % try converting into something numeric
                if(isnumeric(tmp{1}))
                    tmp=cell2mat(tmp');
                end
            catch
                % proceed
            end
            f={field}; % field is a character array, setfield needs a cell array of strings
            dst_data=setfield(dst_data,f{:},tmp);
        end
    end
    
end