A=90;
T=10;
Nperiods=3;
Trisefall=5;
Tsamp=0.1; %% FIXED, Cannot change!!!
Tend=Nperiods*T;
times=0:Tsamp:Tend;
vert=sinetukey(A, T, times, Trisefall);
horiz=sinetukey(A, T, times, Trisefall);
chair=sinetukey(A, T, times, Trisefall);

%horiz=zeros(size(times));
%chair=zeros(size(times));


sr=sr_servo();

sr.write_profile(vert, horiz, chair);
sr.enable();
r=input('press enter to start, q to quit','s');
if ~strcmp(r, 'q')
    sr.start();
    pause(Tend+5);
    %input('press enter to stop')
    sr.stop();
    [pv_vert, pv_horiz, pv_chair]=sr.read_profile_pv();
    [sv_vert, sv_horiz, sv_chair]=sr.read_profile_sv();
    
    times=0:Tsamp:(Tend+5);  % extend plot with 5 seconds
    range=1:length(times);
    vert=zeropad(vert,length(times));
    chair=zeropad(chair,length(times));
    horiz=zeropad(horiz,length(times));
    
    %% Vertical
    figure('Name','Vert Axis');
    title('SV');
    subplot(131)
    plot(times, vert(range),'.');
    title('SV');
    grid('on');
    
    subplot(132);
    plot(times,pv_vert(range),'.');
    title('PV PLC');
    grid('on');
    
    subplot(133);
    plot(times,sv_vert(range),'.');
    title('SV PLC');
    grid('on');
    
    %% Horizontal
    figure('Name','Horiz Axis');
    title('SV');
    subplot(131)
    plot(times, horiz(range),'.');
    title('SV');
    grid('on');
    
    subplot(132);
    plot(times,pv_horiz(range),'.');
    title('PV PLC');
    grid('on');
    
    subplot(133);
    plot(times,sv_horiz(range),'.');
    title('SV PLC');
    grid('on');
    
    % Seat
    figure('Name','Seat Axis');
    title('SV');
    subplot(131)
    plot(times, chair(range),'.');
    title('SV');
    grid('on');
    
    subplot(132);
    plot(times,pv_chair(range),'.');
    title('PV PLC');
    grid('on');
    
    subplot(133);
    plot(times,sv_chair(range),'.');
    title('SV PLC');
    grid('on');
    
end

delete(sr);
