%

classdef biox_rz6_client < biox_abstract_client
    
    properties (Access=protected)
        rz6;
    end

    methods (Access=public)
        function this = biox_rz6_client(number, circuit)
            if number == -1
                this.rz6 = rz6_dummy;
                return;
            end                       

            this.rz6=actxcontrol('RPco.x',[5 5 26 26]);
            if ~this.rz6.ConnectRZ6('GB',number)
               error('rz6.ConnectRZ6');
            end
            
            this.rz6.reset();
            pause(0.05); %RCX file needs time to reset            
            if ~this.rz6.LoadCOF(circuit)
               error('LoadCOF');
            end

            this.rz6.Run();
            pause(0.1); %RCX file needs time to start up                       
            this.write('SYS_SampleRate',this.rz6.GetSFreq());
        end

        function delete(this)
           delete(this.rz6);
        end

        function write(this, tagname, value, offset)
            if nargin < 4
                offset = 0;
            end
            if isscalar(value)
                if ~this.rz6.SetTagVal(tagname,value)
                    error('SetTagVal');
                end
            else
                v=reshape(value,[],1);
                if ~this.rz6.WriteTagVEX(tagname, offset, 'F32', v)
                   error('WriteTagVEX');
                end
            end
        end                       

        function data=read(this, tagname, offset, nWords, datatype, nChannels)
            if nargin < 3
                offset = 0;
            end
            if nargin < 4
                nWords = 1;
            end
            if nargin < 5
                datatype = 'F32';
            end
            if nargin < 6
                nChannels=1;
            end                    
            if nWords > 1
              data=this.rz6.ReadTagVEX(tagname,offset,nWords, datatype, datatype ,nChannels);
            else
              data=this.rz6.GetTagVal(tagname);
            end  
        end
                
        function trigger(this, type)            
            type = lower(type);
            switch type
                case 'soft1' 
                    this.rz6.SoftTrg(1);    
                case 'soft2' 
                    this.rz6.SoftTrg(2);    
                case 'soft3'
                    this.rz6.SoftTrg(3); 
                otherwise
                   error('triggertype can be "soft1", "soft2" or "soft3"');            
            end
        end
        
        function reset_list(this)
            %soft4 resets the tasklist
            this.rz6.SoftTrg(4); 
        end
               
    end
end
