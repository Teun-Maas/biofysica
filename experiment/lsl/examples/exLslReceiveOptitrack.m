function otdata = exLslReceiveOptitrack
    %  info=lsl_resolver('type=''Digital Events @ clockpi'' and name=''Digital Events 0''');
    %  info=lsl_resolver('type=''Pupil Capture @ dcn-eyebrain'' and name=''Pupil Primitive Data - Eye 0''');
    

    %
    % Select Pupil Labs event stream
    %
%    info=lsl_resolver('type=''OptiTrack Mocap @ DCN19'' and name=''Labeled Markers''');
    info=lsl_resolver();

    l=info.list();
    if isempty(l)
        error('no streams found');
    end
  %  l={l};   % FIX THIS IN lsl_resolver.list() !!!
    for i=1:size(l,2)
        fprintf('%i: name: ''%s'' type: ''%s''\n',i,l(i).name, l(i).type);
    end
    
    n=input('enter OptiTrack stream number to acquire: ');
    otstr=lsl_istream(info{n});
    inf = otstr.info();% otstr.info();
fprintf('The stream''s XML meta-data is: \n');
fprintf([inf.as_xml() '\n']);
fprintf(['The manufacturer is: ' inf.desc().child_value('manufacturer') '\n']);
fprintf(['The cap circumference is: ' inf.desc().child('cap').child_value('size') '\n']);
fprintf('The channel labels are as follows:\n');
ch=inf.desc();
ch = inf.desc().child('fields');
ch=ch.child('field');
for k = 1:inf.channel_count()
    fprintf(['  ' ch.child_value('label') '\n']);
    ch = ch.next_sibling();
end


    xmldesc=inf.desc();
    xmldesc.value
    ses=lsl_session();
    ses.add_stream(otstr);
    addlistener(otstr,'DataAvailable',@ot_listener);
    
    input('press enter to start');
    
    ses.start();
    input('press enter to stop');
    
    ses.stop();
    otdata=otstr.read();
    delete(ses);
    delete(otstr);
    delete(info)
    
end

function ot_listener(src, event)
    disp('ot_listener called');
    event
end
