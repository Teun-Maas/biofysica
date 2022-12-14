classdef lsl_istream < handle
    % LSL_ISTREAM - this class represents an lsl_inlet and provides a
    % mechanism to easily collect data by using the LSL_SESION class.
    %
    % See also: LSL_SESSION, LSL_RESOLVER
    %
    % The lsl_session object collects data periodically and also notifies
    % the LSL_ISTREAM objects by posting a DataAvailable event.
    % Listen for these events using addlistener:
    %
    % addlistener(the_lsl_istream,'DataAvailable',@listener_fcn);
    % function listener_fcn(src, event) %#ok<INUSL>
    %     disp('listener_fcn called');
    %     event %#ok<NOPRT>
    % end
    %  
    % The 'event' passed to the listener_fcn will be an lsl_eventdata
    % object, containing all the fields of the lsl_data class.
    %
    % See also: lsl_data, lsl_eventdata, notify, handle, handle.addlistener, event.EventData

    events
        DataAvailable
    end
    
    properties (Access=protected)
        isstring
        inlet
        databuffer
    end
    
    methods
        function this=lsl_istream(streaminfo)
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
        
        function inlet_info=info(this)
            % INFO - returns data provided by the associated inlet's
            % inlet.info() function.
            %
            % See also: lsl_inlet.info
            inlet_info = this.inlet.info();
        end
        
        function result=set_postprocessing(this, processing_flags)
            %SET_POSTPROCESSING - Set timestamp correction postprocessing flags. This should
            % only be done in online scenarios. This option will destroy a
            % stream's original, ground-truth timestamps.
            %
            % See also: lsl_inlet.set_postprocessing
            result=this.inlet.set_postprocessing(processing_flags);
        end

    end
    
    
    methods (Access={?lsl_session})
        
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
                ds{i}=d{1}; %#ok<*AGROW>
                ts{i}=t;
                [d,t]=this.inlet.pull_sample(0);
            end
            timestamps=cell2mat(ts);
            data=ds(1:length(timestamps));
        end
        
        function [data,timestamps]=pull_chunk(this)
            % Wrapper to solve segfault problem mentioned above
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
            % This is called periodically by an lsl_session object
            [data,timestamps]=this.pull_chunk();
            if ~isempty(data)
                tc=this.inlet.time_correction();
                tcindex=length(timestamps);
                evdata=lsl_eventdata(data,timestamps,tc,tcindex);
                this.databuffer.append(evdata);
                notify(this, 'DataAvailable', evdata);
            end
        end
        
        function open_stream(this) %#ok<MANU>
        %    this.inlet.open_stream();
        end

        function close_stream(this) %#ok<MANU>
        %    this.inlet.close_stream();
        end

        function flush_stream(this)
            [~,~]=this.pull_chunk();
            delete(this.databuffer);
            this.databuffer=lsl_databuffer();
        end
    end
    
end
    
%     function timecorrection(this,data)
%         x=data.TCindex;
%         tc=data.TimeCorrection;
%         figure(103);
%         clf;
%         plot(x,tc-tc(1),'o');
%         
%         coeffs=polyfit(x,tc,1);
%         xx=1:numel(ts);
%         xx=indexes;
%         fitted_tc=polyval(coeffs,xx);
%         hold on;
%         plot(xx,fitted_tc-tc(1),'.');
%         hold off;
%         
%         corrected_timestamps=tr+fitted_tc;
%     end
