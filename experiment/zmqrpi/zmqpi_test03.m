function [trigdata, evdata] = zmqpi_test03

  trig=zmqrpi_remote_control('lsldert02.local',5555);

  info=lsl_resolver('type=''Digital Triggers @ lsldert02''');
  l=info.list();
  if isempty(l)
      error('no streams found');
  end
   
  for i=1:numel(l)
     fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type); 
  end
  
  n=1; %input('enter stream number to acquire: ');
  
  trigstr=lsl_istream(info{n});

  info=lsl_resolver('type=''Digital Events @ lsldert02''');
  l=info.list();
  if isempty(l)
      error('no streams found');
  end
   
  for i=1:numel(l)
     fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type); 
  end
  
  n=1 ; %=input('enter stream number to acquire: ');
  
  evstr=lsl_istream(info{n});

  ses=lsl_session();
  ses.add_stream(trigstr);
  addlistener(trigstr,'DataAvailable',@triglistener);

  ses.add_stream(evstr);
  addlistener(evstr,'DataAvailable',@evlistener);

  input('press enter to start');

  ses.start();
  
  for dt = [ 1 2 3 1 2 1 ]
     pause(dt);
     fprintf('playing beep\n');
     trig.beep();
  end
  pause(3);  

  input('press enter to stop');
  delete(trig);
  ses.stop();
  trigdata=trigstr.read();
  evdata=evstr.read();

  delete(ses);
  delete(trigstr);
  delete(evstr);
  delete(info)

end

function triglistener(src, event)
   disp('trigger data received');
   event
end

function evlistener(src, event)
   disp('event data received');
   event
end
