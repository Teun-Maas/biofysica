classdef lsl_ostream < handle

properties (Access=protected)
    outlet

end

methods

    function this=lsl_ostream(info)
        this.outlet=lsl_outlet(info);
    end

    function delete(this)
        delete(this.outlet);
    end

    function write(this,data)
        this.outlet.push_chunk(data);
    end

end
    
end
