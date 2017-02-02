%
leds = pa_leds('fp-web2');
dt=3;
leds.trigger_enable(1);
for i=1:100
    
    disp('voor aan');
    leds.set(2,'red', 1);
    leds.flush;
    disp('voor trigger_wait');
    leds.trigger_wait;
    
    disp('voor uit');
    leds.set(2,'red', 0);
    leds.flush;
    disp('voor trigger_wait');
    leds.trigger_wait;   
end

delete(leds);
