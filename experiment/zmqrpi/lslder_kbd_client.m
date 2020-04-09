classdef lslder_kbd_client < lsldert_zmq_subscriber
   
    methods
        function this=lslder_kbd_client(hostname,port)
            this@lsldert_zmq_subscriber(hostname,port);
            this.socket.subscribe('EV_KEY');
         end
    end  
    
end