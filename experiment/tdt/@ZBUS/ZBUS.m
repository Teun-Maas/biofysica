function [module,err,errstr] = ZBUS(nRacks)
%
% constructor of class zBus
%
% number: number of unit
err         = 0;
module      = actxcontrol('ZBUS.x',[1 1 1 1]);
connect     = module.ConnectZBUS('GB'); % connect to zBus
if ~connect
    err = -1;
    errstr = {'zBus failed to connect'};
else
    for i=1:nRacks
        if (module.HardwareReset(i) ~= 0)
            err = -2;
        end
    end
    if ~err
        for i=1:nRacks
            if (module.FlushIO(i) ~= 0)
                err = -3;
            end
        end
    end
end

