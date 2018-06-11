function data = exLslReceiveData1
%  info=lsl_resolver('type=''Digital Events @ clockpi'' and name=''Digital Events 0''');
%  info=lsl_resolver('type=''Pupil Capture @ dcn-eyebrain'' and name=''Pupil Primitive Data - Eye 0''');
  info=lsl_resolver('type=''Digital Events @ lslder00''');

 
  l=info.list();
  if isempty(l)
      error('no streams found');
  end
   
  for i=1:size(l,1)
     fprintf('%d: name: ''%s'' type: ''%s''\n',i,l{i}.name,l{i}.type); 
  end
  
  n=input('enter stream number to acquire: ');
  
  evstr=lsl_istream(info{n});
  ses=lsl_session();
  ses.add_stream(evstr);
  addlistener(evstr,'DataAvailable',@listener);

  input('press enter to start');

  ses.start();
  input('press enter to stop');
  ses.stop();
  data=evstr.read();
  delete(ses);
  delete(evstr);
  delete(info)

end

function listener(src, event)
   disp('listener called');
   event
end
