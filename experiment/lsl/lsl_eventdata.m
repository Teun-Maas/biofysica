classdef (ConstructOnLoad) lsl_eventdata <  lsl_data & event.EventData
    % lsl_eventdata is used in notify() when sending lsl_data objects.
    
    methods
        function this=lsl_eventdata(data,timestamps,timecorrection,tcindex)
            this@lsl_data(data,timestamps,timecorrection,tcindex);
            this@event.EventData();
        end
    end
    
end