classdef lsl_databuffer < handle
    % lsl_databuffer is used to collect incoming data chunks/objects in
    % lsl_istream
    % the merge_data() methos takes care of combining single
    % chuncks into one contiguous lsl_data object.
    % SEE ALSO: APPEND, READ, CLEAR
    
    % IMPORTANT: when writing methods manipulating this buffer, watch out
    % for data corruption. append() is called from 
    % a timer callback and modifies buf, ibuf and bufsz.
    % When changin these properties from a different thread a semaphore
    % should be used to encapsulate the operations.
    % A semaphore implementation might be found here:
    % https://nl.mathworks.com/matlabcentral/fileexchange/45504-semaphore-posix-and-windows
    
    % FIXME: clear() needs a semaphore.
    
    properties (Access = private)
        buf = {};
        ibuf = 0;
        bufsz = 0;
        readpos = 1;
    end
    
    properties
        chunk = 100;
    end

    methods
        function clear(this)
            % CLEAR - clear the buffer content and reset all housekeeping
            % variables. Currently this is not thread safe. Should be
            % protected by a semaphore. So do not call clear() when an
            % lsl_sesion is collecting data.
            this.buf = {};
            this.ibuf = 0;
            this.bufsz = 0;
            this.chunk = 100;
            this.readpos = 1;
        end
        
        function append(this,data)
            % APPEND - append lsl_data object to buffer
            % append(lsldata);
            if this.ibuf >= this.bufsz
                this.extend(this.chunk);
            end
            this.ibuf=this.ibuf+1;
            this.buf{this.ibuf}=data;
        end
        
        function [data,numread]=read(this,ndata)
            % READ - read merged lsl_data objects from buffer
            % 
            if nargin < 2
                ndata=this.ibuf-this.readpos+1;
            end
            last=this.readpos+ndata-1;

            data=this.merge_data(this.readpos,last);
            this.readpos=last+1;
            numread=ndata;
        end
    end
    
    methods (Access=protected)
        function extend(this,chunk)
            this.bufsz=this.bufsz+chunk;
            this.buf{this.bufsz}=[];
        end
        
        function data=merge_data(this,first, last)
            % MERGE_DATA - merge all lsl_data objects from buf{first} to 
            % buf{last} into a single one
            
            % FIXME/CHECK: this might not work with cf_string type data.
            nbuf=last-first+1;
            if nbuf<1
                % no data or only one sample
                data=[];
                return;
            end
            
            % count total number of samples
            n_alldata=0;
            for i=first:last
               n_alldata=n_alldata+length(this.buf{i}.Timestamps);
            end
            mchan=size(this.buf{first}.Data,1); % each row holds a sample
            
            % are we dealing with cell array data?
            do_cell=iscell(this.buf{first}.Data);
            
            % construct target matrixes
            if do_cell
                merged_data=cell(mchan,n_alldata);
            else
                merged_data=zeros(mchan,n_alldata);
            end
            merged_timestamps=zeros(1,n_alldata); 
            merged_timecorrection=zeros(1,nbuf);
            merged_tcindex=zeros(1,nbuf);
            
            % some running indexes
            nn=1;
            tcbase=0;
            
            % loop over data structures and copy into targets
            for i=first:last
               mrows=size(this.buf{i}.Data,2);
%                if do_cell
%                    merged_data{1:mchan,nn:nn+mrows-1}=this.buf{i}.Data;
%                else
                    merged_data(1:mchan,nn:nn+mrows-1)=this.buf{i}.Data;
%                end
               merged_timestamps(nn:nn+mrows-1)=this.buf{i}.Timestamps;
               merged_timecorrection(i)=this.buf{i}.TimeCorrection;
               % indexes to timecorrection values
               tc=this.buf{i}.TCindex;
               merged_tcindex(i)=tcbase+tc;
               tcbase=tcbase+tc;
               nn=nn+mrows;
            end
            
            data=lsl_data(merged_data,merged_timestamps,...
                merged_timecorrection,merged_tcindex);           
        end
    end
end
