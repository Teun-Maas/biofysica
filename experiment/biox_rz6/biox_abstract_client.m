classdef  biox_abstract_client
    
    
    methods (Abstract)
        write(this, tagname, value, offset)
        data = read(this, tagname, offset, nWords, nChannels)
    end
    
    methods
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
        
        function r=read_acqdata(this, chanlist)
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
        
        function r=read_trialready(this)
            %TODO is dit de juiste variabele?
            r=this.read('STM_Ready');
        end
        
        function r=read_samplerate(this)
            %TODO is dit de juiste variabele?
            r=this.read('SYS_SampleRate');
        end
        
        function write_signalbyte(this, b)
            %TODO is dit de juiste variabele?
            this.write('SGN_Byte',b');
        end
        
        function r=read_inputbyte(this)
            %TODO is dit de juiste variabele?
            r=this.read('INP_Byte');
        end
    end
end
