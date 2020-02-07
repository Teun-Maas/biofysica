function trigdata = lsldert_test01
    
    trig=lsldert_client('lsldert01.local',5555);
    
    info=lsl_resolver('type=''Digital Triggers @ lsldert01''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end
    
    for i=1:numel(l)
        fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
    end
    
    n=1; %input('enter stream number to acquire: ');
    
    trigstr=lsl_istream(info{n});
    
    ses=lsl_session();
    
    cleanupObj=onCleanup(@()cleanupFun);
    
    try
        ses.add_stream(trigstr);
        addlistener(trigstr,'DataAvailable',@triglistener);
                
        input('press enter to start');
        
        ses.start();
        % pause(0.5);
        
        %   for dt = [ 1 2 3 1 2 1 ]
        %      pause(dt);
        %      fprintf('playing beep\n');
        %      trig.beep();
        %   end
        %  pause(3);
        trig.set_digitalin_marker('DRUKKNOP');
        trig.beep(1760, 1);
        trig.digitalout(1);
        doutval=1;
        for jj=1:5
            for ii=1:10
                trig.beep(440,0.5,'PIEP');
                trig.digitalout(doutval,'AU');
                doutval=1-doutval;
                pause(0.02);
            end
            pause(0.27);
        end
        
        input('press enter to stop');
        delete(trig);
        ses.stop();
        trigdata=trigstr.read();
        
    catch
        fprintf('\ncaught exeption\n');
    end
    
    function triglistener(~, event)
        disp('trigger data received');
        if numel(event.Data) ~= 10
            event.Data
        end
    end
        
    function cleanupFun()
        delete(ses);
        delete(trigstr);
        delete(info);
    end
end
