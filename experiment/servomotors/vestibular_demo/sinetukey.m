function [ y ] = sinetukey(A, T, t, Ts)
%SINETUKEY Produce a sine function tapered by a tukey window
%   sinetukey(A, T, t, Ts) produces a sine function
%   y=A*sin(2*pi*t./T) that is windowed using a tukey window
%   rising and falling within Ts seconds

y=A*sin(2*pi*t./T);
w=tukeywin(length(t),Ts/max(t))';
y=y.*w;

end
