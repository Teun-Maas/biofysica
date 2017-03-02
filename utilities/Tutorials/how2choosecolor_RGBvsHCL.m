close all hidden
clear all hidden
clc

%%
n = 64;
subplot(212)
col = pa_statcolor(n,[],[],[],'def',2,'disp',true);
title('HCL')

subplot(211)
RGB = jet(n);
RGB = flipud(RGB);
[m,n] = size(RGB);
RGB		= reshape(RGB,1,m,n);
image(RGB)
axis off;
title('RGB')

%%
print('-depsc','-painter',mfilename);
