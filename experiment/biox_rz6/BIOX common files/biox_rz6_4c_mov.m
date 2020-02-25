classdef biox_rz6_4c_mov < biox_rz6_client

    methods
        function this = biox_rz6_4c_mov(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('BIOX_4C_50kHz_MOV.rcx');
            this@biox_rz6_client(rz6number, f);            
            
        end
        
        function mov_spm_dac(this, dac)
          switch lower(dac) 
            case  'out-a' 
                this.write('MOV_sp0_is_A', 1);
            case  'out-b' 
                this.write('MOV_sp0_is_A', 0);                        
          end;    
        end
        
        %RL: variable 'MOV_Sp_Array' op de RZ6 accepteert alleen een array van lengte 21 
        %RL: Als de lijst een oneven lengte heeft --> Error
        %RL: Als de lijst langer is dan 21 rijen --> Error        
        %RL: Als de lijst korter is dan 21 rijen dan wordt hij aangevuld
        %RL: met dummydata            
        %RL: alist : Matrix van MUX_id's en MUX_channels
        
        function mov_sp_list(this, alist)
            
            listsize = size(alist); 
            if rem(listsize(1),2) == 0   %listsize is even
                error('List size should be odd');
            end    
            if listsize(1) > 21
                error('List size should be <=21'); 
            end
            
            MUX_ids      = alist(:,1);
            MUX_channels = alist(:,2);
            data_in = 16*MUX_ids + MUX_channels;
            
            if listsize(1) < 21
              dummydata = zeros([(21 - listsize(1))/2, 1]);                             
              data_out = [dummydata; data_in; dummydata];
            end
            
            this.write('MOV_Sp_Array', data_out);
        end;
    end
end
