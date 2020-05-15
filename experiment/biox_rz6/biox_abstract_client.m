classdef  biox_abstract_client < handle

    properties (Access=protected)
        % scale factors for acq channels
        % ch5-ch10 are from RA8GA; they need about 1750x in order to volts.
       acq_multipliers = [1 1 1 1 1750 1750 1750 1750 1750 1750 1 1 1 1]; %14 channels    
    end
    
     
    methods (Abstract)
        write(this, tagname, value, offset)
        data = read(this, tagname, offset, nWords, nChannels)   
        trigger(this, type) 
        resetlist(this)
    end
    
    methods
        function write_tasklist(this, tasklist)
            x=tasklist.get();            
            this.write('STM_Matrix',x');  
            this.resetlist(); 
        end
        
        function write_DACoffsets(this, Offsets_A, Offsets_B)
            this.write('OFS_DAC-A1',Offsets_A(1));
            this.write('OFS_DAC-A2',Offsets_A(2));
            this.write('OFS_DAC-A3',Offsets_A(3));
            this.write('OFS_DAC-B1',Offsets_B(1));
            this.write('OFS_DAC-B2',Offsets_B(2));
            this.write('OFS_DAC-B3',Offsets_B(3));
        end
        
        function write_wavdata(this, data, chanlist)             
            if (size(data,2) < 3)
              mydata = transpose(data);
            else
              mydata = data;
            end
            
            if (size(mydata,1) ~= length(chanlist))
               error('Dimensions of WAVdata and chanlist are different'); 
            end   
            
            for i=1:length(chanlist)
                chan=chanlist(i);
                nsamp=size(mydata,2);
                sizetag=sprintf('WAV%d_Size',chan);
                datatag=sprintf('WAV%d_Data',chan); 
                this.write(sizetag,nsamp);
                this.write(datatag,mydata(i,:));            
            end
        end
        
        % RL functie read_acqsize toegevoegd
        function r=read_acqsize(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                readytag=sprintf('ACQ%d_Size',chan);
                r(i)=this.read(readytag);             
            end
        end    
        
        function r=read_acqdata(this, chanlist)
            r=cell(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                sizetag=sprintf('ACQ%d_Size',chan);
                datatag=sprintf('ACQ%d_Data',chan);
                szi=this.read(sizetag);
                multiplier = this.acq_multipliers(chanlist(i));
                r{i}= multiplier * this.read(datatag,0,szi);
            end
        end
        
        function r=read_acqready(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                readytag=sprintf('ACQ%d_Ready',chan); 
                r(i)=this.read(readytag); 
            end
        end
        
        function r=read_trialready(this)
            r=this.read('STM_Ready'); 
        end
        
        function r=read_version(this)
            version=this.read('Version'); 
            r = 'BIOX V3.' + string(version);
        end
        
        function r=read_timer(this)
            r=this.read('Timer'); 
        end
        
        function r=read_samplerate(this)
            r=this.read('SYS_SampleRate');
        end
                
        function write_signalbyte(this, b)
            this.write('SGN_Byte',b'); 
        end
        
        function r=read_inputbyte(this)
            r=this.read('INP_Byte');          
        end
        
        function r=read_inputholdbyte(this)
            r=this.read('INP_HoldByte');         
        end
        
        function r=read_responsetime(this)
            r=this.read('INP_Time')/this.read('SYS_SampleRate');         
        end
        
        function r=read_tasklist(this, tl)  
            r=[];
            for i = 1:tl.nr_of_tasks()
              r(i,:) = this.read('STM_Matrix',7*(i-1), 7,1); %#ok<AGROW>
            end  
        end
        
        function r=read_taskindex(this)
            r = 1 + this.read('STM_CurInd')/7;
        end
                                      
    end
end
