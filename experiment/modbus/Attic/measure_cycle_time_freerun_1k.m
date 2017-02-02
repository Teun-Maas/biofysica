% Time the freerun cycle time for writing the complete led map
% this costs typically 23 ms for a leds.set + leds.flush cycle 
% when the rs232 connection from fpweb2 to plc is running at 119.2 kbps
leds = pa_leds('fp-web2');
disp('now running for 1000 cycles');
tic;
for i=1:500
    leds.set(2,'red', 1);
    leds.flush;
      
    leds.set(2,'red', 0);
    leds.flush;
       
end
toc

delete(leds);
