function lslhid_test02
    
    client=lslder_kbd_client('raspi4',5558);
    
    for i=1:10
       % pause(1);
        s=client.recv_str(0);
        if isempty(s)
            s="nothing";
        end
        s
        fprintf("""%s\n""",s);
    end
    
    delete(client);
end