classdef lsl_data < handle
    properties
        Data
        Timestamps
        TimeCorrection
        TCindex
    end
    
    methods
%         function delete(this)
%             disp('lsl_data.delete() called');
%         end
        
        function this=lsl_data(data,timestamps,timecorrection,tcindex)
            this.Data=data;
            this.Timestamps=timestamps;
            this.TimeCorrection=timecorrection;
            this.TCindex=tcindex;
        end
    end
    
end