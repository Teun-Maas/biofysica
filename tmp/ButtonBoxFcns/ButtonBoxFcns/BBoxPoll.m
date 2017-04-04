function [ButtonVec,ET,AllButtons] = BBoxPoll(Device, MaxWait, MultiWait)
% BBoxPoll -- Poll the response box until one or more buttons is
%   pressed. Device is the Sys3 device name (typically RX6). MaxWait
%   is the maximum time in millisec to wait for the first press. 
%   Set to 0 or Inf to wait forever. MultiWait is the time after the 
%   first press to wait for additional presses. (Optional) MultiWait=0 means
%   only wait for one button. ButtonVec is a 4-element vector with the
%   status of buttons 1 through 4. ET is the elapsed time in millisec
%   to the first press. Return ET= 0 if no press within MaxWait.
%   AllButtons is the decimal equivalent of ButtonVec.
%
% USAGE: [ButtonVec,ET,AllButtons] = BBoxPoll(Device, MaxWait, MultiWait)

KILLTIME = 25;       %    // 25 ms wait while polling

t0 = clock; % starting time

if MaxWait <= 0,
   MaxWait = Inf;
end

if nargin<3,
    MultiWait= 0;
end

ET= 0;
%ButtonVec= BBoxRead(Device);
ButtonVec= [0 0 0 0];

% wait for button press
while ~any(ButtonVec) && ET<(MaxWait/1000),
    MSDelay(KILLTIME);
    ButtonVec= BBoxRead(Device);
    ET= etime(clock,t0);
end

if ~any(ButtonVec), % no activity -- timed out
    ET= 0;
    AllButtons= 0;
    return
end

MSDelay(MultiWait); % wait for more button activity
[ButtonVec,AllButtons]= BBoxRead(Device);
ET= etime(clock,t0);

BBoxClear(Device);  

