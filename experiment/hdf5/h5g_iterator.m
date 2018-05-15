classdef h5g_iterator < h5abstract_iterator
    
    methods
        function this=h5g_iterator(obj)
            this=this@h5abstract_iterator(obj,'h5group','Groups');
        end
        
        function obj=getobj(this,ind)
            obj=h5group(this.obj,this.fields(ind).Name);
        end
    end
end