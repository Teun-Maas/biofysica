classdef biox_rz6_3c < biox_rz6_client
    % biox_rz6_3c Class for communication with a 3 core RZ6 (TDT)    
    
    methods
        function this = biox_rz6_3c(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('biox_3c_50khz.rcx');
            this@biox_rz6_client(rz6number, f);        
        end
    end
end
