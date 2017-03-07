function [r_text_data, r_station_number, r_response_name] = m2c_send(handle, station_nr, command_code, text_data)

   len = length(text_data);
   if len <= 109
      % text_data fits into one frame
      send_single_message_frame(handle, station_nr, command_code, text_data);
   else
      % spread text_data over 2 or more frames
      % first frame text_data must be 109
      send_first_message_frame(handle, station_nr, command_code, text_data(1:108));
      receive_data_request_message_frame(handle);
      i = 109;
      while i+111 < len  % last frame can hold up to 112 chars of text_data 
         j = i+107; % 'next' frames hold 108 chars of text_data
         send_next_message_frame(handle, station_nr, text_data(i:j));
      receive_data_request_message_frame(handle);
         i=i+108;
      end
      send_last_message_frame(handle, station_nr, text_data(i:len));
   end
  % [r_text_data, r_station_number, r_response_name] = receive_response_message(handle);
end

function send_single_message_frame(handle, station_nr, command_code, text_data)
   frame = sprintf('%%%.2d#%.2s%s**', station_nr, command_code, text_data);
   m2c_send_frame(handle, frame);
end

function send_first_message_frame(handle, station_nr, command_code, text_data)
   frame = sprintf('%%%.2d#%.2s%s**&', station_nr, command_code, text_data);
   m2c_send_frame(handle, frame);
end

function send_next_message_frame(handle, station_nr, text_data)
   frame = sprintf('%%%.2d%s**&', station_nr, text_data);
   m2c_send_frame(handle, frame);
end

function send_last_message_frame(handle, station_nr, text_data)
   frame = sprintf('%%%.2d%s**', station_nr, text_data);
   m2c_send_frame(handle, frame);
end

function receive_data_request_message_frame(handle)
   % receive data request message frame '%|NR|BCC|&|CR'
   [ frame, length, ~ ] = m2c_recv_frame(handle);
   if (length ~= 6) || (frame(6) ~= '&')
      error('expected data request message frame, but received "%s"', frame);
   end
end

function [text_data, station_number, response_name] = receive_response_message(handle)
   % receive data request message frame '%|NR|$|RESPNAME|RESPTXT|BCC|CR'
   [ text_data, station_number, response_name ] = m2c_recv(handle);
   if (text_data(4) ~= '$')
      error('expected response message frame, but received "%s"', text_data);
   end
end
