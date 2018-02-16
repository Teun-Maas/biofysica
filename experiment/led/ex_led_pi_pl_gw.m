% Sample program sendig matrixes of LED stimuli to the PLC

%javaclasspath('/Users/gunter/Documents/MATLAB/java/jeromq.jar')
import org.zeromq.ZMQ

n = 4;
s = ledpattern(n);

for i=1:2:n-1
    %s(i).set(0:2:31,'r');
    s(i).set(0:31,'g');

    ir=50;
   % ir=i*50.0/7;
    ig=ir;
   % ig=50-ir;
    s(i).intensity('r', ir);
    s(i).intensity('g', ig);
end

leds = ledcontroller_pi('dcn-led00.local','dcn-led01.local');
%leds = ledcontroller_pi('dcn-led00.local');
leds.print_version;

leds.write(s);
for i=1:n
    input('Press Enter to trigger')
    leds.trigger;
end
if leds.wait
    disp('Done, press enter to mop up');
else
    disp('Been waiting too long..., press enter');
end
pause;
delete(leds);
delete(s);
