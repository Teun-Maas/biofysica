function nwritten = m2c_send_frame(handle, frame)

   global m2c_debuglevel
   if exist('m2c_debuglevel','var') & (m2c_debuglevel > 0)
       fprintf('> "%s"\n', frame);
   end
   pnet(handle, 'printf', '%s\r', frame);

end
