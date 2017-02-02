clear all;
MPD = [];
PPD = [];
latency = [];
blink = [];
baseline = [];
meantrace = [];

names = {
    '02'
    '04'
    '06'
    '08'
    '10'
    '14'
    '16'
    '18'
};
for i = 1:8;
path = ['p',names{i},'_celC.mat'];
load(path)
data = preprocessing(rawdata,50);
meantrace = [meantrace; data{1,2}];

MPD = [MPD data{1,8}];
PPD = [PPD data{1,9}];
latency = [latency data{1,10}];
blink = [blink data{1,5}];
baseline = [baseline data{1,7}];

end
trace_celB = mean(meantrace);
timeB = data{1,1};

figure
hold on
plot (timeB, trace_celB,'linestyle','-','color', 'r')
ylim ([-0.1 0.6]);
xlim ([-1.5 6.0]);
line ([-2 8], [0 0], 'color', [0 0 0], 'linestyle','--')
line ([-1 -1], [-0.5 0.6], 'color', [0 0 0], 'linestyle',':')
line ([0 0], [-0.5 0.6], 'color', [0 0 0], 'linestyle',':')
line ([4.4 4.4], [-0.5 0.6], 'color', [0 0 0], 'linestyle',':')
hold off