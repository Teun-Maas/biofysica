%
leds = ledcontroller('fp-web2');
dt=0.1;
for i=1:16
    leds.set(i,'red', 1);
    leds.flush;
    pause(dt);
    leds.set(i,'grn', 1);
    leds.flush;
    pause(dt)
    
    leds.set(i,'red', 0);
    leds.flush;
    pause(dt);
    leds.set(i,'grn', 0);
    leds.flush;
    pause(dt)
end

delete(leds);
