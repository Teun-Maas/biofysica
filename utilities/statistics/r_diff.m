function p = r_diff(R1,N1,R2,N2)

% P = R_DIFF(R1,N1,R2,N2)
%   two-side significance level of difference between 
%   two correlation coeffients R1 and R2 based on 
%   N1 and N2 data points, respectively. (N1>=10, N2>=10)
%   


%% Fisher's Z-transformation
Z1      = 1/2.*log((1+R1)./(1-R1)); 
Z2      = 1/2.*log((1+R2)./(1-R2));  


num     = abs(Z1-Z2);
den     = sqrt(2) .* sqrt( 1./(N1-3) + 1./(N2-3));

%% Chance
p       = erfc(num./den);



