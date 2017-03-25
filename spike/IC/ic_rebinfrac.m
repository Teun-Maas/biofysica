function y1=ic_rebinfrac(y,fac);

%
%	function y1=rebinfrac(y,fac);
%

a=mod(fac,1);
   
if a==0
   y1=ic_rebinint(y,fac,1);
elseif mod(1/a,1)==0
	fac2=round(1/a);
   fac1=fac*fac2;
   z=ic_rebinint(y,fac2,-1);
   y1=ic_rebinint(z,fac1,1);
else
   y1=y;
   disp(' no easy fraction ')
end

%disp(['new length: ',num2str(length(y1))]);