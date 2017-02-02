function [module, err, errstr] = RZ6(number,circuit)
% [MODULE,ERR, ERRSTR] = RZ6(NUMBER,CIRCUIT)
%
% Initialize RZ6 (constructor of class RZ6)
%
% NUMBER	= number of RP2 as defined by zBusMon
% CIRCUIT	= RPvsdEx circuit filename
% MODULE	= Object
% ERR		= any potential errors
% ERRSTR    = Description of errors

err     = 0;
module  = actxcontrol('RPco.x',[1 1 1 1]); % prog_id, fighandle
connect = module.ConnectRZ6('GB',number); % connect the module
if connect
    module.reset;
    load = module.LoadCOF(circuit); % Load the circuit into module
    if (load)
        module.Run(); % Run the circuit on module
                errstr =   {'No error'};
    else
        err = -2; % Fail to load circuit
        errstr =   {'Failed to connect RZ6'};

    end
else
    err     = -1; % Fail to connect RZ6
    errstr = {'Failure to load ' circuit};
end