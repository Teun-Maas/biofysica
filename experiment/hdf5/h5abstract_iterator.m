classdef h5abstract_iterator < handle
    
    properties (Access=protected)
        obj h5location
        infofield
        classtype
        fields
        i
        n
    end
    
    methods
        function this=h5abstract_iterator(obj,classtype,infofield)
            this.infofield=infofield;
            this.classtype=classtype;
            
            this.obj=obj;
            info=obj.info();
            this.fields=info.(this.infofield);
            this.n=length(this.fields);
            this.i=0;
        end
        
        function result=hasnext(this)
            result = (this.i<this.n);
        end
        
        function rewind(this)
            this.i=0;
        end
        
        function n=numel(this)
            n=this.n;
        end
        
        
        function result=getbyind(this,ind)
            if (ind<1)||(ind>this.n)
                result=[];
            else
                result=this.getobj(ind);
            end
        end
        
        function result=next(this)
            if this.i<this.n
                this.i=this.i+1;
                result=this.getobj(this.i);
            else
                result=[];
            end
        end
        
    end
    
    methods (Abstract)
        obj=getobj(this,ind)
    end
    
end