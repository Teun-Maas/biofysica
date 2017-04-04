function [ButtonVec,AllButtons] = BBoxRead(Device)
% BBoxRead -- Read the buttons on the response box. Device is the
%   Sys3 device (normally RX6). ButtonVec is a 4-element vector with 0 or 
%   1 for buttons 1 through 4. AllButtons is the decimal equivalent of ButtonVec.
%
% USAGE: [ButtonVec,AllButtons]= BBoxRead(Device)

for i=1:4,
    ButtonVec(i)= Device.GetTagVal(['Button',num2str(i)]);
end

ButtonVec= (ButtonVec~=0);  % convert values of the tag to logicals

AllButtons= 0;
AllButtons= sum(bitset(AllButtons,find(ButtonVec))); % convert the vector of logicals to an integer
