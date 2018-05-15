classdef h5a_iterator < h5abstract_iterator
    
    properties (Access=protected)

    end
    
    methods
        function this=h5a_iterator(obj)
            this=this@h5abstract_iterator(obj,'h5attribute','Attributes');
        end
        
        function obj=getobj(this,ind)
            obj=h5attribute(this.obj,this.fields(ind).Name);

%             p=this.fields(ind).Name;
%             obj=feval(this.classtype,this.obj,p);
        end
    end
end