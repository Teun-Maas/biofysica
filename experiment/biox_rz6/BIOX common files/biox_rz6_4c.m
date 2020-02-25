classdef biox_rz6_4c < biox_rz6_client

    methods
        function this = biox_rz6_4c(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('BIOX_4C_50kHz.rcx');
            this@biox_rz6_client(rz6number, f);                     
        end
    end
end
