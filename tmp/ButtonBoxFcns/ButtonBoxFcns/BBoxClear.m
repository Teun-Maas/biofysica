function BBoxClear(Device)
% BBoxClear -- Clear any spurious button presses and wait until all buttons released

KILLTIME = 25;       %    // 25 ms wait while polling

while any(BBoxRead(Device)),
    MSDelay(KILLTIME);
    Device.SoftTrg(2); % relatch the inputs
    drawnow; %flush the event cue, i.e. check for GUI button-presses etc
end
