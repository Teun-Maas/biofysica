function ip=gethostbyname(hostname)
% GETHOSTBYNAME -- returns the IP address of HOSTNAME or an empty char
% array if it is not found

% Günter Windau 20200904
try
    ip = char(py.socket.gethostbyname(hostname));
catch
    ip = '';
end

end