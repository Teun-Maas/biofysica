% RL: This code is for any RZ6 3DSP
% RL: This code works with 'RCX_Uni_3C_50kHz_V3.10.rcx' and above

classdef biox_rz6_3c < biox_rz6_client

    methods
        function this = biox_rz6_3c(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('BIOX_3C_50kHz.rcx');
            this@biox_rz6_client(rz6number, f);        
        end
    end
end
