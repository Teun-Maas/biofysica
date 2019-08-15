classdef rz6_unircx_processor < handle
    properties (Access=protected)
        zBus;
        rz6;
    end

    methods
        function this = rz6_unircx_processor(zBus, rz6)
           this.zBus = zBus;
           this.rz6 = rz6;
           rz6.Run;
        end

        function delete(this)
           delete(this.rz6);
           delete(this.zBus);
        end

        function upload_tasklist(this, tasklist)
            x=tasklist.get();
            this.write('STM_Matrix',x);
        end

        function write(this, tagname, value, offset)
            if nargin < 4
                offset = 0;
            end
            if isscalar(value)
                if ~this.rz6.SetTagVal(tagname,value)
                    error('SetTagVal');
                end
            else
                v=reshape(value,[],1);
                if ~this.rz6.WriteTagVEX(tagname, offset, 'I32', v)
                   error('WriteTagVEX');
                end
            end

        end

        function data=read(this, tagname, nWords, offset, nChannels)
            if nargin < 4
                offset = 0;
            end
            if nargin < 5
                nChannels=1;
            end
            data=this.rz6.ReadTagVEX(tagname,offset,nWords,'I32','I32',nChannels);
        end

    end
end
