%  biox_abstract_client  

classdef  biox_abstract_client < handle

    properties (Access=protected)
        my_version = 32;
        % scale factors for acq channels
        % ch5-ch10 are from RA8GA; they need about 1750x in order to translate to volts.
        acq_multipliers = [1 1 1 1 1750 1750 1750 1750 1750 1750 1]; %11 channels        
    end

    methods (Abstract)
        write(this, tagname, value, offset)
        data = read(this, tagname, offset, nWords, datatype, nChannels)           
        trigger(this, type) 
        reset_list(this)
    end
    
    methods
        function write_tasklist(this, tasklist)
        % WRITE_TASKLIST Adds a tasklist to the clientobject
        %   this.write_tasklist(tl)
        
            x=tasklist.get();            
            this.write('STM_Matrix',x');  
            pause(0.001);
            this.reset_list();
            pause(0.001);
        end
        
        function write_buttonbox_echo(this, varargin)
            p = biox_inputParser;
            p.FunctionName = 'write_buttonbox_echo';
            validOnOffs = { 'On', 'Off' };
            validOnOff = @(x) any(validatestring(x, validOnOffs)); 
            validDuration = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',10});
            
            p.addRequired('OnOff', validOnOff);
            p.addOptional('Duration', 1, validDuration);
            p.parse(varargin{:});
            switch lower(p.Results.OnOff)
                case 'on'
                   onoff = true;
                case 'off'   
                   onoff = false;
                otherwise
                   error('invalid onoff type, this is a bug');   
            end   
            duration = p.Results.Duration;
            this.write('bb_echo',onoff);           
            this.write('bb_echo_dur', duration);
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
                if chan == 11
                    datatype = 'I32';                    
                else
                    datatype = 'F32';                    
                end    
                r{i}= multiplier * this.read(datatag,0,szi, datatype);
                
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
        
        function r=read_wavready(this, chanlist)
            r=zeros(1,length(chanlist));
            for i=1:length(chanlist)
                chan=chanlist(i);
                readytag=sprintf('WAV%d_Ready',chan); 
                r(i)=this.read(readytag); 
            end
        end    
        
        function r=read_trialready(this)
            r=this.read('STM_Ready'); 
        end
        
        function r=read_tasktype(this)
            r=this.read('Task_Type')
        end    
        
        function r=read_rcx_version(this)
            version=this.read('Version'); 
            r = 'BIOX RCX V3.' + string(version);
        end
        
        function r=read_biox_version(this)            
            r = 'BIOX Matlab V3.' + string(this.my_version);
        end
        
        function r=read_timer(this)
            r=this.read('Timer'); 
        end
                      
        function write_signalbyte(this, b)
            this.write('SGN_Byte',b'); 
        end
        
        function r=read_samplerate(this)
            r=this.read('SYS_SampleRate');
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
              r(i,:) = this.read('STM_Matrix',7*(i-1), 7,'F32', 1);
            end  
        end
        
        function r=read_taskindex(this)
            r = 1 + this.read('STM_CurInd')/7;
        end
        
        function r=read_ACQ11(this, varargin)
            
            p = biox_inputParser;
            p.FunctionName = 'read_ACQ11';
            expectedChannel = {'A', 'B', 'C', 'STM', 'A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
                                               
            validChannel = @(x) any(validatestring(x, expectedChannel));                        
            validByte = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<',256});
            
            p.addRequired('Channel', validChannel);            
            p.addOptional('Mask', 255, validByte);
            p.parse(varargin{:});
            
            Channel = lower(p.Results.Channel);    
            mask = uint8(p.Results.Mask);
            
            if strcmp(Channel, 'stm')
                Bit8 = 0;
            end
                       
            switch Channel(1)                
                case 'a'
                    startBit = 0;
                case 'b'
                    startBit = 8;
                case 'c'
                    startBit = 16;
                case 's'  
                    startBit = 24;                    
            otherwise
              error('invalid byte type, this is a bug');
            end
            
            if length(Channel) == 2
                bit8 = str2num(Channel(2));
                mask = 255;
            else
                bit8 = -1;
            end    
            
            if bit8 ~= -1
                bitshiftL = 31 - bit8 - startBit;
                bitshiftR = -31;
            else 
                bitshiftL = 24 - startBit;
                bitshiftR = -24;
            end    
            
            acqdata = this.read_acqdata([11]);
            acqdata11 = uint32(acqdata{1});                     
                           
            for i = 1:length(acqdata11)
                byte = uint8(bitshift(bitshift(acqdata11(i), bitshiftL), bitshiftR)); 
                maskedbyte = bitand(byte, mask); 
                data(i) = maskedbyte;
            end                        
            
            if bit8 == -1
              r = data;
            else
              r = boolean(data);   
            end;  
        end;    

                                      
    end
end
