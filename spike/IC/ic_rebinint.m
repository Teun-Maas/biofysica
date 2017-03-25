function y1=rebinint(y,fac,dir);

% function rebinint(y,fac,dir);
% dir: >0: larger bins; <0: smaller bins
% fac: should be integer;
% 23-1-'01

if nargin<3
   dir=1;
end;

Ny=length(y);
dif=fac-1;
hdif=ceil(0.5*dif);

if dir>0,
	Ny1=floor(Ny/fac);
   y1(1)=mean([y(1:hdif) y(end-hdif+1:end)]);
   for i=2:Ny1,
   	j=fac*(i-1)+1-hdif;
   	y1(i)=mean(y(j:j+dif));
   end;
else
   Ny1=floor(Ny*fac);
	for i=1:Ny,
      j=fac*(i-1)+1;
      y1(j:j+dif)=y(i);
   end;
end;




