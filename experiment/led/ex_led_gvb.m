
% Sample program sendig matrixes of LED stimuli to the PLC
for loc_g =0:127
n=2;
s = ledpattern(n);

loc_g
for i=1:2:n
    s(i).set(loc_g,'g');
%     s(i).set(loc_r,'r');
%     s(i).intensity('r', i);
    s(i).intensity('g', n+1-i);
end

leds = ledcontroller;
leds.write(s);

if leds.wait
	disp('Done, press enter to mop up');
else
    disp('Been waiting too long..., press enter');
end
%pause(2);
tic
delete(leds);
delete(s);
toc
end