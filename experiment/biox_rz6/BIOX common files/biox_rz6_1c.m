classdef biox_rz6_1c < biox_rz6_client
    methods
        function this = biox_rz6_1c(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('BIOX_1C_25kHz.rcx');
            this@biox_rz6_client(rz6number, f);             
        end
    end
end
