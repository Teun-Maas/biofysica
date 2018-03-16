function cfg = sphere_speakerpositionsMinor(cfg,varargin)
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
fname		= which('SphereMinor Measurement.xlsx');% To be checked
N			= xlsread(fname,'Values');
channel		= N(:,1);

%% Transform to Cartesian
% To be checked
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
% [dX,dY,dZ]	= azel2cart(daz,del,R);


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


% sel = sign(del)~=sign(ael) & del<-5 &daz>0;
% [channel(sel) daz(sel) del(sel) aaz(sel) ael(sel)];
