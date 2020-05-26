function lslhid_test02
    
    client=lslder_kbd_client('lsldert01',5558);
    
    for i=1:10
        [key,value]=client.getkey(0);
        fprintf("""%s""\t%d\n",key,value);
    end
    
    delete(client);
end