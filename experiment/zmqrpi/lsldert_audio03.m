client=lsldert_audioclient('raspi6');
for ii=1:10
    %fnii=sprintf('%d.wav',ii);
    client.load(ii,'mariobiondi.mp3', [18 20]*44100);
end

for jj=1:1
    for ii=1:10
        client.play(ii);
        pause(3);
        client.stop(ii);
    end
end

%input('press enter >>>>');

delete(client);
