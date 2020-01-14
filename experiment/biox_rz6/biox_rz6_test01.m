tl=biox_rz6_tasklist;
rz6=biox_rz6_ruurd();
tl.debug(false);

tl.add_task(100,'WaitForTrigge','ZbusB');
tl.add_task(101,'WaitForTrigger','External', 2);
tl.add_task(102,'SoundA','Tone',1000,0,0);
tl.add_task(103,'SoundB','Tone',1000,0,0);
tl.add_task(104,'SoundA','Sweep',1000,3,500);
tl.add_task(105,'SoundA','Noise',1000,50);
tl.add_task(106,'SoundA','Ripple',100,50,100);
tl.add_task(107,'SoundA','Wav',1);
tl.add_task(108,'SoundB','B=A','None');
tl.add_task(109,'SoundB','B=A','Linear',2000,-90);
tl.add_task(110,'SoundB','B=A','Sine',2500,-80);
tl.add_task(1110,'SoundA','Stop');
tl.add_task(112,'Mux', 2, 'reset');
tl.add_task(112,'Mux', 2, 'Set', 1);

tl.add_task(114,'SoundMov','Start',21,2000,-90);
tl.add_task(115,'SoundMov','Stop');
tl.add_task(116,'Daq',15,4,'Start');
tl.add_task(116,'Daq',15,2,'Start');

tl.add_task(117,'Daq',15,4,'Stop');
tl.add_task(118,'SetDio',255);
tl.add_task(119,'TrigOUt',2);
tl.add_task(120,'TrigOUt',2,150);
tl.add_task(121,'reset');
tl.add_task(122,'Ready');
tl.add_task(123,'HoldInpu',4);
% 
rz6.write_tasklist(tl);
input('press enter to finish> ');
delete(tl);
delete(rz6);
