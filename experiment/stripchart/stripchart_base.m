classdef stripchart_base < handle
    
    properties
        ax 
        span
    end
    
    methods
        function this=stripchart_base(haxis)
            this.ax = haxis;
            xlim=get(haxis,'XLim');
            this.span = xlim(2)-xlim(1);
        end
        
    end
    
    methods (Abstract)
        clearpoints(this)
        addpoints(this, Timestamps, Ydata)
        
    end
        
end