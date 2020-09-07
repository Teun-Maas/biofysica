function ip = getip(hostname)
%GETIP Get ip address of a computer by name.
%
% ip = getip returns ip for the computer on which MATLAB is running.
%
% ip = getip('hostname') returns ip for the 'hostname'

%   Copyright (C) Peter Volegov 2002, Albuquerque, NM, USA
%
%   Based on the article contributed by Jeff Lundgren
%   (http://www.codeguru.com/network/local_hostname.shtml) and
%   tcp_udp_ip toolbox by Peter Rydesäter (Peter.Rydesater@mh.se)

% Günter Windau: 20200904 added python workaround if mex file is missing

warning('You need to compile getip.c to an mex file for your platform. Read header of getip.c');
warning('Trying to use python2/3 workaround');

ip=gethostbyname(hostname);