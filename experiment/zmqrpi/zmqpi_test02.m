p=zmqrpi_remote_control('lsldert99.local',5559);

tic;
p.beep();
toc

delete(p);