function Sweep = locsmooth(Sweep,Span,Azimuth,Elevation)
% SWEEP = LOCSMOOTH(SWEEP,SPAN,AZIMUTH,ELEVATION)
%

%% Initialization
n=size(Sweep,2);

%%
S = Sweep;
for i = 1:n
	Disparity	= sqrt((Azimuth - Azimuth(i)).^2+(Elevation - Elevation(i)).^2);
	sel			= Disparity <= Span;
	if any(sel)
		S(:,i) = mean(Sweep(:,sel),2);
	end
end
Sweep = S;

