% close all;
t=tic; 
RZ6_1circuit = which('dataacquisition_v2.rcx');

zBus		= ZBUS(1); % zBus, number of racks
RZ6_1		= RZ6(1,RZ6_1circuit); % Real-time acquisition
nsamples = 5000;
delay = 100; % [ms]
RZ6_1.SetTagVal('acqSamples',nsamples); % amount of DA samples
RZ6_1.SetTagVal('delayAcq',delay);

dur = nsamples/(RZ6_1.GetSFreq/100);
zBus.zBusTrigA(0, 0, 2); % reset, clock start, (0,0,2): trigger entire rack, with a pulse structure, and 2 ms delay(2 ms = minimum).
zBus.zBusTrigB(0, 0, 2); % start event 1/trial onset; trigger zBus 4 = RA16;

% Either of the following blocks works: 
% pause(dur);
 
while RZ6_1.GetTagVal('Active') 
	% do nothing
end
% Using both appears redundant. At least in this example.
a = RZ6_1.ReadTagV('Data_1',0,nsamples)';
b = RZ6_1.ReadTagV('Data_2',0,nsamples)';
c = RZ6_1.ReadTagV('Data_3',0,nsamples)';
h1=figure();
subplot(221)
plot(a)
subplot(222)
plot(b)
subplot(223)
plot(c)
subplot(224)
plot(a,b)

toc(t)