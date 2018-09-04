function nirs = mbll(nirs)
% C = MBLL(OD,WL,DST,DPF)
%
% Modified Beer-Lambert Law
% OD_{\lambda} = \epsilon_{\lambda} \dot c \dot L \dot DPF + OD_{R,\lambda}
% Assume  OD_{R,\lambda} (oxygen-independent optical losses) is constant,
% then:
%
% \Deltac = \DeltaOD_{\lambda}/(\epsilon_{\lambda} \dot c \dot L \dot DPF)
%
% See also MOLAREXTCOEF_HEMO

%% Initialization
% dist		= [hdr.optodedistance hdr.optodedistance]; % Distance between source and detector
% DPF		= hdr.DPF; % differential pathway factor DPF
% wavlen	= hdr.laserwavelength
epsilon	= molarextcoef_hemo; % Molar extinction coefficients
epsilon = round(epsilon);

%%
if ~isfield(nirs,'preprOD')
	OD = nirs.OD;
else
	OD = nirs.preprOD;
end
DPF		= nirs.DPF(1);
dist	= [nirs.distance nirs.distance];
wavlen	= round(nirs.wavelengths);
k = 0;
dC		= NaN(size(OD));
for ii = 1:2:size(OD,2)
	k = k+1;
	d		= dist(k)/10*DPF; % in cm
	
	las1	= wavlen(ii); % source laser 1 wavelength (nm)
	las2	= wavlen(ii+1); % source laser 2 wavelength (nm)
	
	eHbO1	= epsilon(epsilon(:,1) == las1,2)/10^6; % (uM)
	eHbR1	= epsilon(epsilon(:,1) == las1,3)/10^6;
	eHbO2	= epsilon(epsilon(:,1) == las2,2)/10^6;
	eHbR2	= epsilon(epsilon(:,1) == las2,3)/10^6;
	
	dOD1	= OD(:,ii) - median(OD(:,ii)); % change in optical density
	dOD2	= OD(:,ii+1)- median(OD(:,ii+1));  % change in optical density
	
	eM		= [eHbR1*d eHbO1*d; eHbR2*d eHbO2*d];
	dOD		= [dOD1';dOD2'];
	dX		= eM^-1*dOD;
	
	dC(:,ii)	= dX(1,:)'; % deoxy
	dC(:,ii+1)	= dX(2,:)'; % oxy
	
	
end


%% Insert O2HB and HHb labels
label	= nirs.label;
nlabel	= size(label,2);
str		= [];
dclabel = cell(nlabel*2,1);
k		= 0;
for ii = 1:nlabel
	for jj = 1:2
		k = k+1;
		switch jj
			case 1
		str = [label{ii} ' O2Hb'];
			case 2
		str = [label{ii} ' HHb'];
		end
		dclabel{k} = str;
	end
end
nirs.dclabel = dclabel;
%%
nirs.dC = dC;
%%
% keyboard