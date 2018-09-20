function otdata = exLslDisplayOptitrack
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

    global haz hel hrot;
    if ~isempty(haz)
        delete(haz);
    end
    if ~isempty(hel)
        delete(hel);
    end
    if ~isempty(hrot)
        delete(hrot);
    end
    figure(103);
    clf;

    axis([0 10 -180 180]);
    grid('on');
    grid('minor');
    haz = animatedline('LineStyle', 'none', 'Marker', '.', 'MarkerSize', 0.1);
    set(haz,'MaximumNumPoints',119*10,'Color','r');
    hel = animatedline('LineStyle', 'none', 'Marker', '.', 'MarkerSize', 0.1);
    set(hel,'MaximumNumPoints',119*10,'Color','g');
    hrot = animatedline('LineStyle', 'none', 'Marker', '.', 'MarkerSize', 0.1);
    set(hrot,'MaximumNumPoints',119*10,'Color','b');

    ses.start(0.2);
    input('press enter to stop');
    
    ses.stop();
    otdata=otstr.read();
    delete(ses);
    delete(otstr);
    delete(info)
    
end

function ot_listener(src, event)
   % disp('ot_listener called');
    global haz hel hrot;
    persistent t0 span numpoints;
    if isempty(t0)
        t0 = event.Timestamps(1);
        span = 10;
        numpoints=1200;

    end
    X=mod(event.Timestamps-t0,span);

    Xbegin=X(1);
    Xend=X(end);
%     if (Xend <= Xbegin)
%       if numpoints > 0
%          set(haz,'MaximumNumPoints',numpoints);
%          set(hel,'MaximumNumPoints',numpoints);
% 
%          % numpoints=0;
%       end
%     end
%     
    qx = event.Data(11,:);
    qy = event.Data(12,:);
    qz = event.Data(13,:);
    qw = event.Data(14,:);
    [az,el,rot] = quaternion2azel(qx,qy,qz,qw);
%     numpoints=numpoints+numel(az);
%     set(haz,'MaximumNumPoints',numpoints);
%     set(hel,'MaximumNumPoints',numpoints);

    addpoints(haz, X, az);
    addpoints(hel, X, el);
    addpoints(hrot, X, rot);

    drawnow limitrate;
%    nmax=get(h,'MaximumNumPoints')
 %   event
end

