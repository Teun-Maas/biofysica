classdef biox_rz6_4c_ripple < biox_rz6_client

    methods
        function this = biox_rz6_4c_ripple(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end            
            
            f = which('biox_4c_50khz_ripple.rcx');
            this@biox_rz6_client(rz6number, f);                                              
        end        
    end    
end
