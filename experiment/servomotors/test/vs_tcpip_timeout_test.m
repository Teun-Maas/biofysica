
delays=[20 40 60];
n=200;

for ii=delays
    vs=vs_servo;
    times=(0:n-1)/10;
    svax1 = zeros(1,n);
    svax2 = zeros(1,n);
    fprintf('write: ');
    tic;
    vs.write_profile(svax1,svax2);
    toc
    tic;
    fprintf('pausing for %d seconds\n', ii);
    pause(ii);
    toc
    fprintf('before stop\n');
 %   m2c_debug(1);
    tic;
%    vs.keepalive();
    toc
   % m2c_debug(0);
    
    vs.stop();
    fprintf('after stop\n');
    m2c_debug(0);
    delete(vs);
end




