clear all;
close all;

listing = dir('D:\Lidwien\matrixtestPhonak\MatlabOlsa\OlsaMaterial\Sounds_NL');
sndAll = [];
%for ii = [11308 12122 55395 62756 91597 94335 95981 81140 64749 29123 37210  58625 55516  77778 73139 88858 89066 87633 87980 85279 64079 52534 19122 38542 38622 20456];
for ii = 5:300
    filename = ['D:\Lidwien\matrixtestPhonak\MatlabOlsa\OlsaMaterial\Sounds_NL\' listing(ii).name]
   [snd fs] = wavread(filename); 
   snd = circshift(snd(1:fs*1.7,1), round(rand(1)*fs*1.6));
    sndAll = [sndAll snd];
    
    plot(snd); hold all
end

sndAll2 = mean(sndAll');

sndAll2 = [sndAll2 sndAll2 sndAll2 sndAll2 sndAll2];

rmsS = rms(sndAll2);
sndAll3 = sndAll2 / (rmsS/0.1);
rms(sndAll3)

wavwrite(sndAll3, fs, 'dutchspeechnoise')