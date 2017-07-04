function y2 = yalign(x1,y1,x2,y2)


	mn		= max([min(x1) min(x2)]);
	mx		= min([max(x1) max(x2)]);
	xi		= linspace(mn,mx,100);
	y1i		= interp1(x1,y1,xi);
	y2i		= interp1(x2,y2,xi);
	b		= regstats(y1i,y2i,'linear','beta');
	y2		= y2.*b.beta(2)+b.beta(1);
	
	