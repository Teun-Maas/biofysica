%
leds = pa_leds('fp-web2');
dt=0.02;
leds.set(2,'red', 1);
leds.set(3,'grn', 1);
leds.flush;

for i=0:100
    leds.intensity('red', i)
    leds.intensity('grn', i)
end
%leds.clear_all
delete(leds);
