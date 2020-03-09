function result=compare_timestamps(data)
   n = numel(data);
   ts = cell(1,n);
   
   for ii=1:n
       ts{ii} = lsl_correct_lsl_timestamps(data{ii});
   end
   offset=ts{1}(1);
   for ii=1:n
       ts{ii} = (ts{ii}-offset);
   end
   
   m=numel(ts{1});
   result=zeros(m,n);
   for ii=1:n
       result(:,ii) = ts{ii};
   end   
   
end