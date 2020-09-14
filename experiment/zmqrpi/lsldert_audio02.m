client=lsldert_audioclient('raspi5.local');
for ii=1:10
    fnii=sprintf('%d.wav',ii);
    client.load(ii,fnii);
end

for jj=1:2
    for ii=1:10
        client.play(ii);
        pause(1.5);
        client.stop(ii);
    end
end

%input('press enter >>>>');

delete(client);
