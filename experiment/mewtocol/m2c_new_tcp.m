function handle = m2c_new_tcp(hostname, port)

   handle = pnet('tcpconnect', hostname, port);
   pnet(handle,'setreadtimeout', 5)


end

