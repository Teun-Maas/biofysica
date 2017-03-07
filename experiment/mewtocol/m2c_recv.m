function [ text_data, station_number, command_response_code ] = m2c_recv(handle)
   %
   %
   %
   [ frame, length, more ] = m2c_recv_frame(handle);
   station_number = frame(2:3);
   bcc = frame(length-1:length);
   command_response_symbol = frame(4);
   command_response_code = frame(5:6);
   if more
       text_data = frame(7:length-3);
   else
       text_data = frame(7:length-2);
   end
   
   if command_response_symbol == '!'
       m2c_error(text_data, station_number, command_response_code);
   end

   while more
      request_next_frame(handle, station_number);
      [ frame, length, more ] = m2c_recv_frame(handle);
      if more
          more_data = frame(4:length-3);
      else
          more_data = frame(4:length-2);
      end    
      text_data = [ text_data more_data ];
   end
end

function request_next_frame(handle, station_number)
   msg = [ '%' station_number '**&' ];
   m2c_send_frame(handle, msg);
end
