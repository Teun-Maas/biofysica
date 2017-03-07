function [ frame, length, more ] = m2c_recv_frame(handle)
   %
   % receive one MEWTOCOL-COM frame
   %
   idx=[];
   while isempty(idx)
      peek = pnet(handle, 'read', 1024, 'char', 'view', 'noblock');
      CR = char(13);
      idx = find(peek==CR, 1);
      length = idx-1;
   end
   frame = pnet(handle, 'read', idx);

   if ((frame(length) == '&') & (length ~= 6))   % length 6 is a data request message frame '%|NR|BCC|&|CR'
      more = 1;
   else
      more=0;
   end
   frame(idx) = [];
   global m2c_debuglevel
   if exist('m2c_debuglevel','var') & (m2c_debuglevel > 0)
      fprintf('< "%s"\n', frame);
   end
end

