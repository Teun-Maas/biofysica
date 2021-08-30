classdef lsl_session < handle
    % LSL_SESSION - session class that handles incoming data from
    % lsl_istream objects in the background.
    %
    % See also: ADD_STREAM, START, STOP
    %
    % LSL_SESSION uses a timer to periodically check attached lsl_istream
    % objects for available data. On errors when the LSL_SESSION objects
    % are not deleted, the timer keeps on running. When creating new
    % sessions and attached streams this may cause unexpected results.
    %
    % Neat solution: perform exception handling in your code that cleans up
    % after errors.
    % Workaround: manually remove all timers running in the background
    % using the command 'delete(timerfindall)'
    % Todo: upon creation check if an LSL_SESSION is still running and act
    % upon this.
    
    properties
        tmr
        lsl_lib
        streams={}
        nstreams=0
    end
    
    methods (Access=public)
        function this=lsl_session()
            % Initialize the lsl_session
            this.lsl_lib=lsl_loadlib();
            this.setup_timer();
        end
        
        function delete(this)
             this.stop();
             delete(this.tmr);
        end        
        
        function add_stream(this,stream)
            % add an lsl_istream object to this session
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
            % start() starts data collection from the attached lsl_istream
            % objects.
            % start(interval_seconds) starts data collection with the
            % specified interval.
            if nargin>1
                this.tmr.period=period;
            end
            this.flush_streams();
            this.open_streams();
            start(this.tmr);
            pause(0.1); 
        end
        
        function stop(this)
            stop(this.tmr);
            this.collect();  % get remaining data
            this.close_streams();
        end
           
        function open_streams(this)
            i=1;
            while i<=this.nstreams
                this.streams{i}.open_stream();
                i=i+1;
            end
        end

        function close_streams(this)
            i=1;
            while i<=this.nstreams
                this.streams{i}.close_stream();
                i=i+1;
            end
        end

        function flush_streams(this)
            i=1;
            while i<=this.nstreams
                this.streams{i}.flush_stream();
                i=i+1;
            end
        end

    end
    
end
