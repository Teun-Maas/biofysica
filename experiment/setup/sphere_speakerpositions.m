function cfg = sphere_speakerpositions(cfg,varargin)
% CFG = SPHERE_SPEAKERPOSITIONS
%
% Convert sperical coordinates to double polar coordinates for the SPHERE lab
% Data is obtained from 'Sphere Measurement.xlsx' with values recorded by
% Sebastian Ausili.
%
% See also SPHERELOOKUP

if nargin<1
	cfg = [];
end

%% Optional arguments
dspFlag       = keyval('display',varargin);
if isempty(dspFlag)
	dspFlag	= false;
end
if dspFlag
	close all
end
%% File
fname		= which('Sphere Measurement.xlsx');
N			= xlsread(fname,'Values');
channel		= N(:,1);

%% Transform to Cartesian
% actual spherical azimuth, elevation
az			= -N(:,4)+90;
el			= N(:,5);
az			= pa_deg2rad(az);
el			= pa_deg2rad(el);
R			= 1;
[X,Y,Z]		= sph2cart(az,el,R);
[X,Y,Z]		= pitch(X,Y,Z,90);

% desired double-polar azimuth and elevation
daz			= N(:,2);
del			= N(:,3);
[dX,dY,dZ]	= azel2cart(daz,del,R);

if dspFlag
	figure
	plot3(X,Y,Z,'.')
	hold on
	plot3(dX,dY,dZ,'o')
	sel		= daz==0;
	[~,I]	= sortrows([daz(sel) del(sel)]);
	x		= [dX(sel) dY(sel) dZ(sel)];
	x		= x(I,:);
	plot3(x(:,1),x(:,2),x(:,3),'k-');
	[~,I]	= sortrows([daz(sel) del(sel)]);
	x		= [X(sel) Y(sel) Z(sel)];
	x		= x(I,:);
	plot3(x(:,1),x(:,2),x(:,3),'r-');
	axis square
	grid on
	xlabel('X');
	ylabel('Y');
	zlabel('Z');
	
	sel		= del==0;
	[~,I]	= sortrows([daz(sel) del(sel)]);
	x		= [dX(sel) dY(sel) dZ(sel)];
	x		= x(I,:);
	plot3(x(:,1),x(:,2),x(:,3),'k-');
	
	[~,I]	= sortrows([daz(sel) del(sel)]);
	x		= [X(sel) Y(sel) Z(sel)];
	x		= x(I,:);
	plot3(x(:,1),x(:,2),x(:,3),'r-');
	title('Cartesian');
end
%% Transform to double-polar coordinate system
% actual double-polar azimuth, elevation
[X,Y,Z]		= yaw(X,Y,Z,-90);
[aazel]		= xyz2azel(X,Y,Z);
aaz			= aazel(:,1);
ael			= aazel(:,2);

sel			= daz>90;
aaz(sel)	= 90+(90-aaz(sel));

sel			= daz<-90;
aaz(sel)	= -90+(-90-aaz(sel));

sel			= del>90;
ael(sel)	= 90+(90-ael(sel));


cfg.lookup	= [aaz ael channel];


if dspFlag
	figure
	subplot(131)
% 	plot(daz,del,'r.');
	hold on
	plot(aaz,ael,'o','MarkerFaceColor','w');
	axis([-130 130 -130 130]);
	axis square;
	
	subplot(232)
	plot(daz,aaz,'.');
	axis square;
	axis([-130 130 -130 130]);
	title('Azimuth');
	ylabel('Actual Position (deg)');
	xlabel('Desired Position (deg)');
	
	subplot(233)
	plot(del,ael,'.');
	axis square;
	axis([-130 130 -130 130]);
	title('Elevation');
	ylabel('Actual Position (deg)');
	xlabel('Desired Position (deg)');
	
	subplot(235)
	plot(daz,aaz-daz,'.');
	axis square;
	axis([-130 130 -60 60]);
	title('Azimuth');
	ylabel('Error (deg)');
	xlabel('Desired Position (deg)');
	
	subplot(236)
	plot(del,ael-del,'.');
	axis square;
	axis([-130 130 -60 60]);
	title('Elevation');
	ylabel('Error (deg)');
	xlabel('Desired Position (deg)');
end

sel = sign(del)~=sign(ael) & del<-5 &daz>0;
[channel(sel) daz(sel) del(sel) aaz(sel) ael(sel)];
