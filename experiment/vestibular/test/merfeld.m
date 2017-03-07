function [ y ] = merfeld(A, T, t, Ts1, Ts2)
%MERFELD Generate a modulated cosine function
% As of 2017-03-03 this function is deprecated because of possible
% discontinuities occuring around Ts1 and Tmax-Ts2
% Please use SINETRPZ which generates a trapezoid tapered sine wave
% 
    if nargin  < 4
        Ts1 = 0;
    end
    if nargin < 5
        Ts2 = 0;
    end
    
    y=zeros(size(t));
    y(t<0) = 0;
    t2 = t(t>=0 & t<Ts1);
    y(t>=0 & t<Ts1) = merfeld_start(A, T, t2, Ts1);
    
    t3 = t(t>=Ts1);
    y(t>=Ts1) = -A*cos(2*pi*t3./T);
    
    tm=max(t);
    t4=t(t>tm-Ts2 & t<=tm);
    y(t>tm-Ts2 & t<=tm) = -merfeld_start(A, T, t4-tm, Ts2);
end

function y = merfeld_start(A, T, t, Ts)
    y = (A/Ts)*((T/(2*pi))*sin(2*pi*t./T)-t.*cos(2*pi*t./T));
end
    
