function trigdata = lsldert_pubsubtest01
    
    trig = lsldert_pubclient;
    
    info=lsl_resolver('type=''Digital Triggers @ lsldert00''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end   
    fprintf('%d: name: ''%s'' type: ''%s''\n',1,l(1).name,l(1).type);
    trigstr1=lsl_istream(info{1});  
    
    info=lsl_resolver('type=''Digital Triggers @ lsldert03''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end
    fprintf('%d: name: ''%s'' type: ''%s''\n',1,l(1).name,l(1).type);
    trigstr2=lsl_istream(info{1});  

    ses=lsl_session();
    
    cleanupObj=onCleanup(@()cleanupFun);
    
%    try
        ses.add_stream(trigstr1);
        ses.add_stream(trigstr2);

        addlistener(trigstr1,'DataAvailable',@triglistener);
        addlistener(trigstr2,'DataAvailable',@triglistener);
                
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
        for jj=1:2
            for ii=1:10
                trig.digitalout(doutval,'AU');
                trig.beep(440,0.05,'PIEP');

                doutval=1-doutval;
                pause(0.1);
            end
            %pause(0.27);
        end
        
        input('press enter to stop');
        delete(trig);
        ses.stop();
        trigdata=trigstr1.read();
%         evdata=evstr.read();
        
%    catch ME
%        fprintf('\ncaught exeption\n');
%        ME
%   end
    
    function triglistener(src, event)
        disp('trigger data received');
        info=src.info();
        info.name
        info.type
        info.source_id
    
        event.Data
    end
        
    function cleanupFun()
        delete(ses);
        delete(trigstr1);
        delete(trigstr2);
        delete(info);
    end
end
