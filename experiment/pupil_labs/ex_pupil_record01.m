%javaclasspath('/usr/share/java/jeromq.jar');

%rc = pupil_remote_control('dcn-eyebrain.local');
%rc = pupil_remote_control('pupil-desktop.local');
rc = pupil_remote_control('localhost');

r=rc.time_sync(0.0);

ntimes=1;
rt=zeros(1,ntimes);
for i=1:ntimes
    tstart=tic;
    timestamp=rc.get_time_stamp();
    roundtrip=toc(tstart)*1000;
    rt(i)=roundtrip;
    data=rc.record;
  %  fprintf(1,"pupil time stamp = %18.15f roundtrip time=%4.1f ms\n",timestamp, roundtrip);
end

r=rc.stop_recording;

delete(rc);
