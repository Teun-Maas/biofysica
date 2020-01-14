classdef biox_rz6_client < biox_abstract_client
    
    properties (Access=protected)
        zBus;
        rz6;
    end

    methods
        function this = biox_rz6_client(number, circuit)
            if number == -1
                this.zBus = zbus_dummy;
                this.rz6 = rz6_dummy;
                return;
            end

            this.zBus = actxcontrol('ZBUS.x',[1 1 1 1]);
            if ~this.zBus.ConnectZBUS('GB')
               error('zBus.ConnectZBUS');
            end

            this.rz6=actxcontrol('RPco.x',[5 5 26 26]);
            if ~this.rz6.ConnectRZ6('GB',number)
               error('rz6.ConnectRZ6');
            end
            
            this.rz6.reset();

            if ~this.rz6.LoadCOF(circuit)
               error('LoadCOF');
            end

            this.rz6.Run();
        end

        function delete(this)
           delete(this.rz6);
           delete(this.zBus);
        end

        function write(this, tagname, value, offset)
            if nargin < 4
                offset = 0;
            end
%TODO check of scalars ook met WriteTagVEX goed worden geschreven
            %if isscalar(value)
            %    if ~this.rz6.SetTagVal(tagname,value)
            %        error('SetTagVal');
            %    end
            %else
                v=reshape(value,[],1);
                if ~this.rz6.WriteTagVEX(tagname, offset, 'I32', v)
                   error('WriteTagVEX');
                end
            %end
        end

        function data=read(this, tagname, offset, nWords, nChannels)
            if nargin < 3
                offset = 0;
            end
            if nargin < 4
                nWords = 1;
            end
            if nargin < 5
                nChannels=1;
            end
            data=this.rz6.ReadTagVEX(tagname,offset,nWords,'I32','I32',nChannels);
        end

    end
end
