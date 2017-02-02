function [module,err,errstr] = RA16(number,circuit)
% [MODULE,ERROR] = RA16(NUMBER,CIRCUIT)
%
% constructor of class RA16
%
err  = 0;
module = actxcontrol('RPco.x',[1 1 1 1]);
connect = module.ConnectRA16('GB',number); % Connect to RA16 via Gigabit
if connect
	module.reset;
	load = module.LoadCOF(circuit); % load circuit
	if load
		module.Run(); % run circuit
		module.SetTagVal('Run',0); % MW: Does this still exist
		%or
		% RA16.ClearCOF;
	else
		err = -2;
		errstr = {'RA16 failed to load'};
	end
else
	err = -1;
	errstr = {'RA16 failed to connect'};
end
