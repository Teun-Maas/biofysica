classdef lslder_kbd_client < lsldert_zmq_subscriber
   
    methods
        function this=lslder_kbd_client(hostname,port)
            this@lsldert_zmq_subscriber(hostname,port);
            this.socket.subscribe('EV_KEY');
        end
         
        function [key, value] = getkey(this,dontwait)
            if dontwait
                str = this.recv_str(1);
            else
                str = [];
                while isempty(str)
                    str = this.recv_str(0);                    
                end
            end
            if isempty(str)
                key='';
                value=NaN;
                return
            end
            ss=strsplit(string(str));
            key=ss{2};
            value=str2double(ss{3});
        end
    end  
    
end