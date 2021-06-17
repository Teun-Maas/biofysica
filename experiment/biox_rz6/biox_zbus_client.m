%  biox_zbus_client

classdef biox_zbus_client < handle
    
    properties (Access=protected)
        zBus;
    end

    methods
        function this = biox_zbus_client()                     
            this.zBus = actxcontrol('ZBUS.x',[1 1 1 1]);
            if ~this.zBus.ConnectZBUS('GB')
               error('zBus.ConnectZBUS');
            end
        end

        function delete(this)
           delete(this.zBus);
        end
        
        function trigger(this, aTrigger)            
            switch lower(aTrigger)
                case 'a'
                    this.zBus.zBusTrigA(0, 0, 5); % all rack, pulse, 5 ms delay
                case 'b'     
                    this.zBus.zBusTrigB(0, 0, 5); % all rack, pulse, 5 ms delay            
            end        
        end
    end
end    
