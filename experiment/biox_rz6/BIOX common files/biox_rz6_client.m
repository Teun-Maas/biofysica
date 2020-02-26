classdef biox_rz6_client < biox_abstract_client
    
    properties (Access=protected)
        rz6;
    end

    methods
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
            if nWords > 1
              data=this.rz6.ReadTagVEX(tagname,offset,nWords,'F32','F32',nChannels);
            else
              data=this.rz6.GetTagVal(tagname);
            end  
        end
        
        %RL: method toegevoegd.
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
                   error('unknown triggertype');            
            end
        end
        
        %RL: method startlist toegevoegd voor resetten van tasklist.
        %RL: deze functie wordt gebruikt in de functie <write_tasklist> in biox_abstract_client
        %RL: maar heb ik bewust ook public gemaakt.  
        function resetlist(this)
            this.rz6.SoftTrg(4); %soft4 resets the tasklist
        end
    end
end
