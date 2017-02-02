%
% constructor of class RA8GA
%
% number: number of unit
function [module error] = RA8GA(number)
    error  = 0;
    module = actxcontrol('RPco.x',[1 1 1 1]);
    connect = module.ConnectRA8GA('GB',number);
    if (connect == 0) 
        error = -1;
    else
        module.reset;
    end
end
