function [veltrace,smvtrace,acctrace,smatrace]          = getvel(htrace,vtrace,Fsample,sd)
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% [VEL, SVEL] = GETVEL(HOR,VER,SMOOTHFACTOR,FSAMPLE)
%
% Obtain radial velocity from horizontal HOR and vertical VER traces.
%
% See also GSMOOTH, PA_SACDET
%
% MarcW 2007

Rx                                      = htrace;
Ry                                      = vtrace;
R                                       = NaN*Rx;
veltrace                                = R;
acctrace                                = R;
smvtrace                                = R;
smatrace                                = R;
for i                                   = 1:size(htrace,2)
	Rx(:,i)                             = gradient(Rx(:,i),1);
	Ry(:,i)                             = gradient(Ry(:,i),1);
	R(:,i)                              = hypot(Rx(:,i),Ry(:,i));
	R(:,i)                              = cumsum(R(:,i));
	
	veltrace(:,i)                       = gradient(R(:,i),1./Fsample);
	smvtrace(:,i)                       = pa_gsmooth(veltrace(:,i),Fsample,sd);
	
	acctrace(:,i)                       = gradient(smvtrace(:,i),1./Fsample);
	smatrace(:,i)                       = pa_gsmooth(acctrace(:,i),Fsample,sd);
end
