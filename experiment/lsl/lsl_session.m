classdef lsl_session < handle
    properties
        tmr
        lsl_lib
        streams={}
        nstreams=0
    end
    
    methods (Access=public)
        function this=lsl_session()
            this.lsl_lib=lsl_loadlib();
            this.setup_timer();
        end
        
        function delete(this)
             this.stop();
             delete(this.tmr);
        end        
        
        function add_stream(this,stream)
            this.nstreams=this.nstreams+1;
            this.streams{this.nstreams}=stream;
        end
        
        function collect(this)
            % let all streams collect data
            % when data is available the stream will trigger a notification
            % event 'DataAvailable' to its listeners.
            i=1;
            while i<=this.nstreams
                this.streams{i}.collect();
                i=i+1;
            end
        end
        
        function setup_timer(this)
            this.tmr = timer('ExecutionMode','fixedRate', ... % Run continuously
                'TimerFcn',@lsl_timer_callback); % Run callback function at every timer event
            this.tmr.period=1.0;
            
            function lsl_timer_callback(~,~)
                this.collect();
            end
            
        end
        
        function start(this,period)  % you can actually stop and restart the nested timer function at any time
            if nargin>1
                this.tmr.period=period;
            end
            this.flush_streams();
            start(this.tmr);
        end
        
        function stop(this)
            stop(this.tmr);
            this.collect();  % get remaining data
        end
           
        function flush_streams(this)
% TODO: make non blocking
%             i=1;
%             while i<=this.nstreams
%                 this.streams{i}.flush();
%                 i=i+1;
%             end
        end
    end
    
end