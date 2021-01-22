classdef biox_rz6_bb_direct < biox_rz6_client

    methods
        function this = biox_rz6_bb_direct(rz6number)
            if nargin < 1
                rz6number=-1;  % use rz6_dummy for debugging
            end
            f = which('biox_3c_BB_direct.rcx');
            this@biox_rz6_client(rz6number, f);
            this.bb_version = 'BB';
        end
    end
end
