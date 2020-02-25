% RL: This code is for RZ6 ID #3960 (3DSP) alone
% RL: This code works with 'RCX_Uni_3C_50kHz_V3.10.rcx' alone

classdef biox_rz6_4c_ripple < biox_rz6_client

    methods
        function this = biox_rz6_4c_ripple(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end            
            
            f = which('BIOX_4C_50kHz_Ripple_V3.16.rcx');
            this@biox_rz6_client(rz6number, f);                                              
        end        
    end
    
    
    
end
