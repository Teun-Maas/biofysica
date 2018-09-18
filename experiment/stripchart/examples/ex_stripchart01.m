% ex_stripchart01

fsamp=120;

scfig=1034;
figure(scfig);
clf(scfig);
ax=axes;%(scfig);
set(ax,'XLim',[0 10], 'YLim', [-2 2]);
sc = stripchart(ax,1200);


chunksz=20;
freq=0.55;
timestamps=(0:chunksz-1)/fsamp;

t0=0;

for k=1:1000
    y=sin(2*pi*freq*timestamps);
    timestamps=timestamps+chunksz/fsamp;
    sc.addpoints(timestamps,y);
    pause(chunksz/fsamp);
end
