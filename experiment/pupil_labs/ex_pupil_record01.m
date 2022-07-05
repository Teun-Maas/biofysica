%javaclasspath('/usr/share/java/jeromq.jar');

%rc = pupil_remote_control('dcn-eyebrain.local');
%rc = pupil_remote_control('pupil-desktop.local');
rc = pupil_remote_control('localhost');

%% DON'T DO THIS WHEN USING LSL 
%% r=rc.time_sync(0.0);

disp('press a key to start recording');
pause;
r=rc.start_recording;

disp('press a key to stop recording');
pause;
r=rc.stop_recording;

delete(rc);
