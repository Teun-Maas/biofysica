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
    end
end
