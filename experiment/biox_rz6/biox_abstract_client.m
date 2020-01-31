classdef  biox_abstract_client

    properties (Access=protected)
       acq_multipliers = [1 1 1 1 1 1 1 1 1 1] ;
    end
    
    
    methods (Abstract)
        write(this, tagname, value, offset)
        data = read(this, tagname, offset, nWords, nChannels)   
        trigger(this, type) 
        startlist(this)
    end
    
    methods
        function write_tasklist(this, tasklist)
            x=tasklist.get();
            this.write('STM_Matrix',x'); %RL: OK   
            %RL startlist() toegevoegd omdat tasklist soms niet wil starten
            this.startlist(); 
        end
        
        %RL: functie write_DACoffsets toegevoegd 
        function write_DACoffsets(this, Offsets_A, Offsets_B)
            this.write('OFS_DAC-A1',Offsets_A(1)); %RL: Naam OK
            this.write('OFS_DAC-A2',Offsets_A(2)); %RL: Naam OK
            this.write('OFS_DAC-A3',Offsets_A(3)); %RL: Naam OK
            this.write('OFS_DAC-B1',Offsets_B(1)); %RL: Naam OK
            this.write('OFS_DAC-B2',Offsets_B(2)); %RL: Naam OK
            this.write('OFS_DAC-B3',Offsets_B(3)); %RL: Naam OK
        end
        
        function write_wavdata(this, data, chanlist)
            for i=1:length(chanlist)
                chan=chanlist(i);
                nsamp=size(data,2);
                sizetag=sprintf('WAV%d_Size',chan); %RL: Naam aangepast 
                datatag=sprintf('WAV%d_Data',chan); %RL: Naam aangepast 
                this.write(sizetag,nsamp);
 
                this.writeF(datatag,data(i,:)); %RL 'chan' verander naar 'i'
                %this.writeF(datatag,data(chan,:)); 
            end
        end
        
        % RL functie read_acqsize toegevoegd
        function r=read_acqsize(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                readytag=sprintf('ACQ%d_Size',chan); %RL: Naam OK
                r(i)=this.read(readytag);             
            end
        end    
        
        function r=read_acqdata(this, chanlist)
            r=cell(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                sizetag=sprintf('ACQ%d_Size',chan); %Naam aangepast RL
                datatag=sprintf('ACQ%d_Data',chan); %Naam aangepast RL
                szi=this.read(sizetag);
                multiplier = this.acq_multipliers(chanlist(i));
                r{i}= multiplier * this.read(datatag,0,szi);
            end
        end
        
        function r=read_acqready(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                readytag=sprintf('ACQ%d_Ready',chan); %RL: Naam aangepast
                r(i)=this.read(readytag); 
            end
        end
        
        function r=read_trialready(this)
            r=this.read('STM_Ready'); %RL: Naam OK 
        end
        
        function r=read_samplerate(this)
            r=this.read('SYS_SampleRate'); %RL: Naam OK 
        end
                
        function write_signalbyte(this, b)
            this.write('SGN_Byte',b'); %RL: Naam OK 
        end
        
        %RL: method toegevoegd.
        function r=read_inputbyte(this)
            r=this.read('INP_Byte'); %RL: Naam OK           
        end
        
        %RL: method toegevoegd.
        function r=read_inputholdbyte(this)
            r=this.read('INP_HoldByte'); %RL: Naam OK           
        end
        
        %RL: method toegevoegd.
        function r=read_responsetime(this)
            r=this.read('INP_Time')/this.read('SYS_SampleRate'); %RL: Naam OK           
        end
                                      
    end
end
