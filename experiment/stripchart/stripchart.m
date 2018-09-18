classdef stripchart < stripchart_base
    % simple strip chart recorder based on animatedline
    % this one assumes a fixed sample rate for incoming data
    properties (Access=protected)
        hline 
        numpoints
        t0
    end
    
    methods
        function this = stripchart(haxis, numpoints)
           this = this@stripchart_base(haxis);
           this.numpoints=numpoints;
           this.t0=[];
           this.hline = animatedline(this.ax, ...
               'LineStyle', 'none',...
               'Marker', '.',...
               'MarkerSize', 0.5,...
               'MaximumNumPoints', this.numpoints...
               );
        end
        
        function clearpoints(this)
            clearpoints(this.hline);
        end
        
        function addpoints(this, Timestamps, Ydata)
            if isempty(this.t0)
                this.t0=Timestamps(1);
            end
            Xdata=mod(Timestamps-this.t0,this.span);

            addpoints(this.hline, Xdata, Ydata);
        end
        
    end
    
    
end