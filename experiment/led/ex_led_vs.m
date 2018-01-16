% Sample program sendig matrixes of LED stimuli to the PLC
% We need jeromq here! Add the full path to jeromq.jar to you javaclasspath.txt
% Adapted for vestibular chair setup with five LEDs

import org.zeromq.ZMQ

% NOTICE: only use hostnames of modules in your own setup!!!
host1='dcn-led02.local';

n = 8;
s = ledpattern(n);

for i=1:5
    s(i).set(i-1,'r');   
    s(i).set(i-1,'g');
    
    ir=50;
    ig=ir;
    % ig=50-ir;   % this will generate warnings
    s(i).intensity('r', ir);
    s(i).intensity('g', ig);
end

leds = ledcontroller_pi(host1);          %one LED box
%leds = ledcontroller_pi(host1, host2);  %two LED boxes in parallel

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
