function MUX(RP,Device,Channel)
% MUXset(RP,DEVICE,CHANNEL)
%
% Set the PM2Relay multiplexers 
%   RP: handle to RP2 object
%   DEVICE: PM2Relay device number 1, 2, 3 or 4. 
%   CHANNEL:  channel number 1 to 16. 
%
%	To disable all channels on a particular device, pass Channel=0 (by
%	default)

if nargin<3
	Channel = 0;
end

DeviceTable = [0 16 32 48];
RP.SetTagVal('DeviceSelect',DeviceTable(Device));  % select the device
if Channel, % activate a channel
    RP.SetTagVal('ChanSelect',Channel-1); % select the channel
    RP.SetTagVal('SetReset',64);    % make it active
else % inactivate all channels on this device
    RP.SetTagVal('SetReset',128);
end
RP.SoftTrg(1);
