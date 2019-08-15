classdef rz6_unircx_1core_50kHz < rz6_unircx_processor

    methods
        function this = rz6_unircx_1core_50kHz
            zBus = actxcontrol('ZBUS.x',[1 1 1 1]);
            if ~zBus.ConnectZBUS('GB')
               error('zBus.ConnectZBUS');
            end

            rz6=actxcontrol('RPco.x',[5 5 26 26]);
            if ~rz6.ConnectRZ6('GB',1)
               error('rz6.ConnectRZ6');
            end

            if ~rz6.LoadCOF('RCX_Uni_1C_50kHz_V3.00.rcx')
               error('LoadCOF');
            end

            this@rz6_unircx_processor(zBus,rz6);
            
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
