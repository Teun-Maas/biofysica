% Sample program sendig matrixes of LED stimuli to the PLC
% We need jeromq here! Add the full path to jeromq.jar to you javaclasspath.txt

%import org.zeromq.ZMQ

% NOTICE: only use hostnames of modules in your own setup!!!

n = 8;
s = ledpattern(n);

for i=1:2:n-1
    s(i).set(0,'r');
    s(i).set(1,'g');
    ir=i*50.0/7;
    ig=ir;
    % ig=50-ir;   % this will generate warnings
    s(i).intensity('r', ir);
    s(i).intensity('g', ig);
end
leds = ledcontroller_pi('dcn-led09','dcn-led10','dcn-led06','dcn-led07');

leds.print_version;
% leds.write(s);
% for i=1:n
%     input('Press Enter to trigger')
%     leds.trigger;
% end
% if leds.wait
%     disp('Done, press enter to mop up');
% else
%     disp('Been waiting too long..., press enter');
% end
pause;
delete(leds);
delete(s);
