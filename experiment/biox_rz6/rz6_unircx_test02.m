%rz6=rz6_unircx_client(1,'RCX_Uni_1C_50kHz_V2.13.rcx');
tl=rz6_unircx_tasklist();
rz6=rz6_unircx_client(-1);

tl.debug(true);

tl.add_task(100,'WaitForTrigger','ZbusB');
%tl.add_task(101,'WaitForTrigger','External', 2);
%tl.add_task(1000,'SoundA','Tone',1000,0,0,0);
%tl.add_task(1200,'SoundB','Tone',1000,0,0,0); 
%tl.add_task(104,'SoundA','Sweep',1000,3,500,0);
%tl.add_task(105,'SoundA','Noise',10000,50,0);
tl.add_task(106,'SoundA','Ripple',100,50,100);
%tl.add_task(107,'SoundA','Wav',6,0);
%tl.add_task(108,'SoundA','B=A','None',0);
%tl.add_task(109,'SoundA','B=A','Linear',2000,-90,0);
%tl.add_task(110,'SoundA','B=A','Sine',2500,-80,0);
tl.add_task(10000,'SoundA','Stop');
%tl.add_task(1500,'SoundB','Stop');
%tl.add_task(112,'Mux', 8);
%tl.add_task(113,'Signaling', 9);
%tl.add_task(114,'SoundMov','Start',32,2000,-90);
%tl.add_task(115,'SoundMov','Stop');
%%tl.add_task(116,'Daq','Start');
%tl.add_task(116,'Daq','Start',15);
%tl.add_task(116,'Daq','Start',15,2);
%
%tl.add_task(117,'Daq','Stop');
%tl.add_task(118,'SetDio',255);
%tl.add_task(119,'TrigOUt',2);
%tl.add_task(120,'TrigOUt',2,150);
%tl.add_task(121,'Reset');
%tl.add_task(122,'Ready');
%tl.add_task(123,'HoldInp','Start',4);
%tl.add_task(124,'HoldInp','Stop');
%
rz6.write_tasklist(tl);
input('press enter to stop> ');
delete(tl);
delete(rz6);
