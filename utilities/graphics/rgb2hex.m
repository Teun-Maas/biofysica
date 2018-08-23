function H = rgb2hex(RGB)
% H = RGB2HEX(RGB)
%
% Convert RGB color values to hexadecimal colour values H
%
% See also PA_STATCOLOR

% 2014 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com


n = size(RGB,1);
H = cell(n,1);
for row = 1:n
	rgb = RGB(row,:);
	rgb = rgb*255;
	hex1 = dec2hex(floor(rgb/16));
	hex2 = dec2hex(round(mod(rgb,16)*16));
	hex = ['#',hex1(1),hex2(1),hex1(2),hex2(2),hex1(3),hex2(3)];
	H{row} = hex;
end