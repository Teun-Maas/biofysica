% Sample program sendig matrixes of LED stimuli to the PLC

n = 8;
s = ledpattern(n);

for i=1:2:n
    s(i).set(0:2:127,'r');
    s(i).set(1:2:127,'g');
    s(i).intensity('r', i);
    s(i).intensity('g', n+1-i);
end

leds = ledcontroller;
leds.write(s);
if leds.wait
    disp('Done, press enter to mop up');
else
    disp('Been waiting too long..., press enter');
end
pause;
delete(leds);
delete(s);
