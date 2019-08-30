classdef rz6_unircx_1core_50kHz < rz6_unircx_client

    methods
        function this = rz6_unircx_1core_50kHz

            this@rz6_unircx_client(1,'RCX_Uni_1C_50kHz_V2.13.rcx');
            
            MUX0=1*16;
            MUX1=2*16;
            SPKMOV=zeros(1,21);
            SPKMOV(1:2:21)=MUX0;
            SPKMOV(2:2:20)=MUX1;
            SPKMOV=SPKMOV+[0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10];
            this.write('MOV_SSp_Array', SPKMOV);
            this.write('MOV_Sp0_is_A',1);
        end
    end
end
