function result=compare_timestamps(data)
   n = numel(data);
   ts = cell(1,n);
   
   for ii=1:n
       ts{ii} = lsl_correct_lsl_timestamps(data{ii});
     %  ts{ii} = data{ii}.Timestamps;
   end

   
   m=numel(ts{1});
   result=zeros(m,n);
   for ii=1:n
       result(:,ii) = ts{ii};
   end   
   
   offset=min(result(1,:));
   result=result-offset;
end

function dingetjes
    ts=compare_timestamps(data);
    
    dts=1000*(ts - ts(:,1));
    plot(dts(:,1))
    plot(ts(:,1),'.')
    plot(data{1}.TCindex,data{1}.TimeCorrection,'.')
    plot(data{1}.Timestamps,'.')
    tc=lsl_estimate_timecorrection(data{2})
end