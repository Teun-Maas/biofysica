function y2 = yalign(x1,y1,x2,y2)
% Y2 = YALIGN(X1,Y1,X2,Y2)
%
% Align curves 1 and 2, by shifting curve 2 in the y-dimension
%
% Alignment is done via interpolation over overlapping x regions and
% regressing out the offset.
%
% See also INTERP1, REGSTATS

mn		= max([min(x1) min(x2)]);
mx		= min([max(x1) max(x2)]);
xi		= linspace(mn,mx,100);
y1i		= interp1(x1,y1,xi);
y2i		= interp1(x2,y2,xi);
b		= regstats(y1i,y2i,'linear','beta');
y2		= y2.*b.beta(2)+b.beta(1);

