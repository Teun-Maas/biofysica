function [module,err,errstr] = PA5(number)
% PA5(NUMBER)
%
% constructor of class PA5
%
% List of PA5 methods
%
% GetError - Retrieve an error message
% str = PA5.GetError;
% if length(str) > 0
%   disp str
% end
%
% GetAtten - Returns current level of attenuation
% atten = PA5.GetAtten;
%
% SetAtten - Sets the attenuation
% PA5.SetAtten(atten)
%
% SetUser - Set User Attenuation mode
% PA5.SetUser(mode, value)
% --> MODE      --> Value
% base      = 1     5   (is used when several speakers
%                        vary in signal intensity)
% step      = 2     0.0 -> 120.0 dB
% reference = 3     0.0 -> 120.0 dB (reference dB level
% update    = 4
%                   0-dynamic(4,0)atten is changed as the dial is turned
%                   1-manual (4,1)updates after pressing select button
% absmin    = 5     20.0  Absolute minimal attenuation
err         = 0;
module      = actxcontrol('PA5.x',[1 1 1 1]);
connect     = module.ConnectPA5('GB',number);
if ~connect
	err = -1;
    errstr = {'PA5 failed to connect'};
else
	module.reset; % sets level to standard 0.0 dB attenuation;

end



