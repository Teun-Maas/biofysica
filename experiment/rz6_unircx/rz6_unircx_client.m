classdef rz6_unircx_client < handle
    properties (Access=protected)
        zBus;
        rz6;
    end

    methods
        function this = rz6_unircx_client(number, circuit)
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

        function write_tasklist(this, tasklist)
            x=tasklist.get();
            this.write('STM_Matrix',x');
        end

        function write_wavdata(this, data, chanlist)
            for i=1:length(chanlist)
               chan=chanlist(i);
               nsamp=size(data,2);
               sizetag=sprintf('BufferSize%d',chan);
               datatag=sprintf('WavData%d',chan);
               this.write(sizetag,nsamp);
               this.write(datatag,data(chan,:));
            end
        end


        function read_acqdata(this, chanlist)
            r=cell(1,length(chanlist));
            for i=1:length(chanlist)
               chan=chanlist(i);
               sizetag=sprintf('Acq%d_Size',chan);
               datatag=sprintf('Acq%d_Data',chan);
               szi=this.read(sizetag);
               r{i}=this.read(datatag,0,szi);
            end
        end

        function r=read_acqready(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
               chan=chanlist(i);
               readytag=sprintf('Acq%d_Ready',chan);
               r(i)=this.read(readytag);
            end
        end

        function write(this, tagname, value, offset)
            if nargin < 4
                offset = 0;
            end
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
