function trigdata = lsldert_pubsubtest04
    
    proxy=lsldert_pubclient('lsldert00');
    ses=lsl_session();
 
    %lslclients = { 'lsldert00', 'raspi4', 'raspi5', 'raspi6'};
    lslclients = { 'lsldert05', 'lsldert00'};

    numclients = numel(lslclients);
    
    for ii=1:numclients
        info=lsl_resolver(sprintf('type=''Digital Triggers @ %s''', lslclients{ii}));
        l=info.list();
        if isempty(l)
            error('no streams found');
        end
        
        for i=1:numel(l)
            fprintf('%d: name: ''%s'' type: ''%s''\n',i,l(i).name,l(i).type);
        end
                
        trigstr(ii)=lsl_istream(info{1});
        ses.add_stream(trigstr(ii));

    end
    
    %cleanupObj=onCleanup(@()cleanupFun);
    
    try
        addlistener(trigstr(1),'DataAvailable',@triglistener); %trigstr(1) is the zmq proxy
        
        input('press enter to start');
        
        ses.start();
        % pause(0.5);
        
        %   for dt = [ 1 2 3 1 2 1 ]
        %      pause(dt);
        %      fprintf('playing beep\n');
        %      trig.beep();
        %   end
        pause(3);
        proxy.set_digitalin_marker('DRUKKNOP');
    %    trig.beep(1760, 1);
        proxy.digitalout(1);
        doutval=1;
        for jj=1:50
            %for ii=1:100
                %  trig.beep(440,0.5,'PIEP');
                proxy.pulseIR(3,5)
             %   proxy.digitalout(doutval,'AU');
                doutval=1-doutval;
                pause(5.5);
            %end
            %pause(0.27);
        end
        
        input('press enter to stop');
        ses.stop();
        
        for ii=1:numclients
            trigdata{ii}=trigstr(ii).read();
        end
        
        delete(ses);
        delete(trigstr);
        delete(info);
        delete(proxy);

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
%         delete(ses);
%         delete(trigstr);
%         delete(info);
    end
end
