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

        end

        function write(this, tagname, value)

        end

        function read(this, tagname, value)

        end

    end
end
