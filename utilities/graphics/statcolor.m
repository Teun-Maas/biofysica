function RGB = statcolor(ncol,statmap,palette,Par,varargin)
% C = STATCOLOR(NCOL,STATMAP,PALETTE)
%
% Choose a color palette for statistical graphs.
%
% Statistical map:
%	- Qualitative
%		for categories, nominal data. Default palettes include:
%			# Dynamic (d)	- 1
%			# Harmonic (h)	- 2 
%			# Cold (c)		- 3
%			# Warm (w)		- 4
%
%	- Sequential
%		for numerical information, for metric & ordinal data, when low
%		values are uninteresting and high values interesting. Default
%		palettes include:
%			# Luminance, l	- 5
%			# LumChrm, lc	- 6
%			# LumChrmH, lch - 7
%
%	- Diverging				- 8
%		for numerical information, for metric & ordinal data, when negative
%		(low) values and positive (high) values are interesting and a
%		neutral value (0) is insignificant. This map has no palette. 
%							
%
%
% Example usage:
% >> ColMap = pa_statcolor(64,'sequential','l',260);
% >> colormap(ColMap)
%
% 	 ColMap = statcolor(64,[],[],[],'def',8)
% 		
% This function is based on:
% Zeileis, Hornik, Murrel. Escaping RGBland: Selecting colors for
% statistical graphics. Computational Statistics and Data Analysis 53
% (2009) 3259–3270
% http://www.sciencedirect.com/science/article/pii/S0167947308005549 
%
% Note: has not been tested
%
% See also:
% http://research.stowers-institute.org/efg/R/Color/Chart/
% http://geography.uoregon.edu/datagraphics/color_scales.htm
% http://hclcolor.com/

% As per point 5 from
% http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_updates.html, a
% citation:
% Brewer, Cynthia A., 2013. http://www.ColorBrewer.org, accessed 29-07-2013

% 2013 Marc van Wanrooij
% e: marcvanwanrooij@neural-code.com

%% Check
if nargin<1
	ncol = 2^8;
	close all
end
if nargin<2
	% 	statmap = 'qualitative';
	statmap = 'sequential';
end
if nargin<3
	palette = 'luminancechroma';
end
if nargin<4
	% 	Par = [Lmin Lmax Cmin Cmax Hmin Hmax];
	Par = [100 0 100 20];
end
dispFlag = keyval('disp',varargin,false);
def = keyval('def',varargin);

%% Default values
switch def
	case 1
		statmap = 'qualitative';
		palette = 'dynamic';
		Par		= [];
	case 2
		statmap = 'qualitative';
		palette = 'harmonic';
		Par		= [];
	case 3
		statmap = 'qualitative';
		palette = 'cold';
		Par		= [];
	case 4
		statmap = 'qualitative';
		palette = 'warm';
		Par		= [];
	case 5
		statmap = 'sequential';
		palette = 'luminance';
		Par		= [];
	case 6
		statmap = 'sequential';
		palette = 'luminancechroma';
		Par		= [0 100 100 20]; % [Lmin Lmax Cmax H]
	case 7
		statmap = 'sequential';
		palette = 'luminancechromahue';
		Par = [0 100 100 100 30 90]; % [Lmin Lmax Cmin Cmax H1 H2] % Heat
	case 8
		statmap = 'diverging';
		palette = [];
		Par = [10 100 100 260 30]; % [Lmin Lmax Cmax H1 H2] 'Blue-White-Red'
		% 		Par = [40 100 90 140 320]; % [Lmin Lmax Cmax H1 H2] 'Green-White-Purple'
	case 9
		statmap = 'sequential';
		palette = 'luminancechromahue';
		Par = [70 70 70 70 0 360]; % [Lmin Lmax Cmin Cmax H1 H2] %Rainbow
	case 10
		statmap = 'divergingskew';
		palette = [];
		Par = [00 100 100 260 30]; % [Lmin Lmax Cmax H1 H2] 'Blue-White-Red'
end
statmap = lower(statmap);
palette = lower(palette);

switch statmap
	case 'qualitative'
		switch palette
			case {'d','dynamic'}
				H		= linspace(30,300,ncol);
			case {'h','harmonic'}
				H		= linspace(60,240,ncol);
			case {'c','cold'}
				H		= linspace(270,150,ncol);
			case {'w','warm'}
				H		= linspace(90,-30,ncol);
		end
		C		= repmat(70,1,ncol);
		L		= repmat(70,1,ncol);
		LCH		= [L;C;H]';
		
		
	case 'sequential'
		switch palette
			case {'l','luminance'}
				l		= linspace(0,2,ncol);
				H		= repmat(300,1,ncol);
				C		= zeros(1,ncol);
				L		= 90-l*30;
				LCH		= [L;C;H]';
				
			case {'lc','luminancechroma'}
				l		= linspace(0,1,ncol);
				p		= 1;
				fl		= l.^p;
				Cmax	= Par(3);
				Lmax	= Par(2);
				Lmin	= Par(1);
				H		= Par(4);
				H		= repmat(H,1,ncol);
				C		= zeros(1,ncol)+fl*Cmax;
				L		= Lmax-fl*(Lmax-Lmin);
				LCH		= [L;C;H]';
				
			case {'lch','luminancechromahue'}
				l		= linspace(0,1,ncol);
				p		= 1;
				fl		= l.^p;
				Cmax	= Par(4);
				Cmin	= Par(3);
				Lmax	= Par(2);
				Lmin	= Par(1);
				H1		= Par(5);
				H2		= Par(6);
				H		= H2-l*(H2-H1);
				C		= Cmax-fl*(Cmax-Cmin);
				L		= Lmax-fl*(Lmax-Lmin);
				LCH		= [L;C;H]';
				
		end
	case 'diverging'
		nhalf	= ceil(ncol/2);
		l		= linspace(0,1,nhalf);
		p		= 1;
		fl		= l.^p;
		Cmax	= Par(3);
		Lmax	= Par(2);
		Lmin	= Par(1);
		H1		= Par(4);
		H2		= Par(5);
		
		H1		= repmat(H1,1,nhalf);
		H2		= repmat(H2,1,nhalf);
		C		= zeros(1,ncol/2)+fl*Cmax;
		L		= Lmax-fl*(Lmax-Lmin);
		LCH1	= flipud([L;C;H1]');
		LCH2	= [L;C;H2]';
		LCH		= [LCH1; LCH2];

	case 'divergingskew'
		nthird	= ceil(ncol/3);
		l1		= linspace(0,1,nthird);
		l2		= linspace(0,1,ncol-nthird);
		p		= 1;
		fl1		= l1.^p;
		fl2		= l2.^p;
		Cmax	= Par(3);
		Lmax	= Par(2);
		Lmin	= Par(1);
		H1		= Par(4);
		H2		= Par(5);
		
		H1		= repmat(H1,1,nthird);
		H2		= repmat(H2,1,ncol-nthird);
		C1		= zeros(1,nthird)+fl1*Cmax;
		C2		= zeros(1,ncol-nthird)+fl2*Cmax;
		L1		= Lmax-fl1*(Lmax-Lmin);
		L2		= Lmax-fl2*(Lmax-Lmin);
		LCH1	= flipud([L1;C1;H1]');
		LCH2	= [L2;C2;H2]';
		LCH		= [LCH1; LCH2];
				
end
RGB		= LCH2RGB(LCH);

if dispFlag
plotcolmap(RGB)
end

function plotcolmap(RGB)
% PLOTCOLMAP(RGB)
% Plot the colors on a scale
% ncol	= size(RGB,1);
% for ii	= 1:ncol
% 	plot(ii,1,'ks','MarkerFaceColor',RGB(ii,:),'MarkerSize',15,'MarkerEdgeColor','k');
% 	hold on
% end
% xlim([0 ncol+1]);
% ylim([0 2])
% axis off

[m,n] = size(RGB);
RGB		= reshape(RGB,1,m,n);
image(RGB)
axis off;