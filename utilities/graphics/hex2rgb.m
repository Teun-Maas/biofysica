function rgb = hex2rgb(hex)
% H = RGB2HEX(RGB)
%
% Convert RGB color values to hexadecimal colour values H
%
% See also PA_STATCOLOR

% 2014 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

if nargin<1
hex = '#A3586D';
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

% hex1(1)=H(2)
% hex2(1)=H(3)
% hex1(2)=H(4)
% hex2(2)=H(5)
% hex1(3)=H(6)
% hex2(3)= H(7)

% dec1 = 16*hex2dec(hex1)
% dec2 = hex2dec(hex2)

rgb = reshape( sscanf(hex,'%2x') ,3,[]).'/255;

% n = size(RGB,1);
% H = cell(n,1);
% for row = 1:n
% 	rgb = RGB(row,:);
% 	rgb = rgb*255;
% 	hex1 = dec2hex(floor(rgb/16));
% 	hex2 = dec2hex(round(mod(rgb,16)*16));
% 	hex = ['#',hex1(1),hex2(1),hex1(2),hex2(2),hex1(3),hex2(3)];
% 	H{row} = hex;
% end
% #A3586D