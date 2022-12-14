function [evdata, pldata] = exLslReceiveDataSphere1
    %  info=lsl_resolver('type=''Digital Events @ clockpi'' and name=''Digital Events 0''');
    %  info=lsl_resolver('type=''Pupil Capture @ dcn-eyebrain'' and name=''Pupil Primitive Data - Eye 0''');
    
    %
    % Select digital events stream
    %
    info=lsl_resolver('type=''Digital Events @ lslder04'' and name=''Digital Events 1''');
    l=info.infos();
    if isempty(l)
        error('no streams found');
    end
    
    for i=1:size(l,1)
        fprintf('%d: name: ''%s'' type: ''%s''\n',i,l{i}.name,l{i}.type);
    end
    
    n=input('enter digital event stream number to acquire: ');
    evstr=lsl_istream(info{n});
    
    %
    % Select Pupil Labs event stream
    %
    info=lsl_resolver('type=''Pupil Capture @ pupil-hpdesktop'' and name=''Pupil Primitive Data - Eye 0''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end
    
    for i=1:size(l,1)
        fprintf('%d: name: ''%s'' type: ''%s''\n',i,l{i}.name,l{i}.type);
    end
    
    n=input('enter digital event stream number to acquire: ');
    plstr=lsl_istream(info{n});
    
    
    
    ses=lsl_session();
    ses.add_stream(evstr);
    ses.add_stream(plstr);
    addlistener(evstr,'DataAvailable',@ev_listener);
    addlistener(plstr,'DataAvailable',@pl_listener);
    
    input('press enter to start');
    
    ses.start();
    input('press enter to stop');
   % pause(30);
    ses.stop();
    evdata=evstr.read();
    pldata=plstr.read();
    delete(ses);
    delete(evstr);
    delete(plstr);
    delete(info)
    
end

function ev_listener(src, event)
    disp('ev_listener called');
    % event
end

function pl_listener(src, event)
    disp('pl_listener called');
    % event
end
