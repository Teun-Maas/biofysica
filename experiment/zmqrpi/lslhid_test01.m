function HIDdata=lslhid_test01
info=lsl_resolver('type=''HID Events @ raspi4''');
    l=info.list();
    if isempty(l)
        error('no streams found');
    end
    fprintf('%d: name: ''%s'' type: ''%s''\n',1,l(1).name,l(1).type);
    HIDstr=lsl_istream(info{1});  

    ses=lsl_session();
            ses.add_stream(HIDstr);

        addlistener(HIDstr,'DataAvailable',@triglistener);
                
        input('press enter to start');
        
        ses.start();

        
        input('press enter to stop');
        ses.stop();
        HIDdata=HIDstr.read();
    
    function triglistener(src, event)
        disp('trigger data received');
        info=src.info();
        info.name
        info.type
        info.source_id
    
        event.Data
    end
        
end