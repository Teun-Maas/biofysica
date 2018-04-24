classdef lsl_stream < handle
    events
        DataAvailable
    end
    
    properties % (Access=protected)
        isstring
        inlet
        databuffer
    end
    
    methods
        function this=lsl_stream(streaminfo)
            this.inlet=lsl_inlet(streaminfo);
            this.inlet.time_correction();  % initialize network time correction protocol
            this.isstring=strcmp(streaminfo.channel_format,'cf_string');
            this.databuffer=lsl_databuffer();
        end
        
        function delete(this)
            delete(this.inlet);
            delete(this.databuffer);
        end
        
        function data=read(this)
            % READ - read all captured data in buffer
            % notice that the processing functions in lsl_databuffer are
            % not thread safe, so do not call read() while
            % the lsl_session is running! calls lsl_session.stop() first.
            data=this.databuffer.read();
        end
        
        function [data,timestamps]=pull_sample(this,timeout)
            [data,timestamps]=this.inlet.pull_sample(timeout);
        end
        
        function [data,timestamps]=pull_string_chunk(this)
            % Workaround. Strings apparently can only be pulled one at a time.
            % pull_chunk() segfaults on them...
            bufsz=0;
            ds={};
            ts={};
            i=0;
            [d,t]=this.inlet.pull_sample(0);
            while ~isempty(d)
                i=i+1;
                if i>bufsz
                    % expand buffer by 1000 cells
                    bufsz=bufsz+1000;
                    ds{bufsz}=[];
                    ts{bufsz}=[];
                end
                ds{i}=d;
                ts{i}=t;
                [d,t]=this.inlet.pull_sample(0);
            end
            timestamps=cell2mat(ts);
            data=ds(1:length(timestamps));
        end
        
        function [data,timestamps]=pull_chunk(this)
            if this.isstring
                [data,timestamps]=this.pull_string_chunk();
            else
                [data,timestamps]=this.inlet.pull_chunk();
            end
        end
        
        function collect(this)
            % try to collect data
            % if data available pack into lsl_eventdata
            % and notify
            [data,timestamps]=this.pull_chunk();
            if ~isempty(data)
                tc=this.inlet.time_correction();
                tcindex=length(timestamps);
                evdata=lsl_eventdata(data,timestamps,tc,tcindex);
                this.databuffer.append(evdata);
                notify(this, 'DataAvailable', evdata);
            end
        end
    end
    
    
end