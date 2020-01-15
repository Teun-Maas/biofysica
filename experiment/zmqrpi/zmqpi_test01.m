p=zmqrpi_remote_control('lsldert99.local',5555);

tic;
p.beep();
toc

delete(p);