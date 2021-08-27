function [pldata, metadata] = exLslReceivePLPy_v2
    %  info=lsl_resolver('type=''Digital Events @ clockpi'' and name=''Digital Events 0''');
    %  info=lsl_resolver('type=''Pupil Capture @ dcn-eyebrain'' and name=''Pupil Primitive Data - Eye 0''');
    

    %
    % Select Pupil Labs event stream
    %
    info=lsl_resolver('type=''Pupil Gaze @ LAPTOP-CBLJOAUJ''');

    l=info.list();
    if isempty(l)
        error('no streams found');
    end
    
    for i=1:size(l,1)
        fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
    end
    
    n=input('enter pupil labs stream number to acquire: ');
    plstr=lsl_istream(info{n});
    
    metadata = lsl_metadata_gaze(plstr);
    %fprintf([metadata.as_xml() '\n']);
    p=metadata.as_struct(); %#ok<NASGU>
    
    ses=lsl_session();
    ses.add_stream(plstr);
    addlistener(plstr,'DataAvailable',@pl_listener);
    
    input('press enter to start');
    
    ses.start();
    input('press enter to stop');
    
    ses.stop();
    pldata=plstr.read();
    delete(ses);
    delete(plstr);
    delete(info);    
end

function pl_listener(src, event) %#ok<INUSL>
    disp('pl_listener called');
    event %#ok<NOPRT>
end
