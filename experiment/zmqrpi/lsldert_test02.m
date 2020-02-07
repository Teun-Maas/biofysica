function trigdata = lsldert_test01
    
    trigc=cell(1);
    trigc{1}=lsldert_client('lsldert01.local',5555);
    trigc{2}=lsldert_client('lsldert04.local',5555);

    trig = lsldert_cluster;
    %trig.add_client(trigc{1});
    %trig.add_client(trigc{2});
    trig.add_client(trigc);
    
    info=lsl_resolver('type=''Digital Triggers @ lsldert01''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end   
    fprintf('%d: name: ''%s'' type: ''%s''\n',1,l(1).name,l(1).type);
    trigstr1=lsl_istream(info{1});  
    
    info=lsl_resolver('type=''Digital Triggers @ lsldert04''');
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
