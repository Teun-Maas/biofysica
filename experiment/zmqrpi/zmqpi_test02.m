p=zmqrpi_remote_control('lsldert02.local',5555);

tic;
p.beep();
p.beep(880, 1.1);
toc

delete(p);