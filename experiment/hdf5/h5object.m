classdef h5object < h5location
        
    methods
        function this = h5object(parent_location,name)
            this=this@h5location(parent_location,name);
        end
        
        %%% Attribute open/create  stuff goes here
        function attr=getattr(this,name)
            attr=h5attribute(this,attr_name);
        end
        
        function attr=createattr(this,attr_name,value)
            attr=h5attribute(this,attr_name);
            if nargin >= 3
                attr.write(value);
            end
        end
        
        function members=info(this)
            members=h5info(this.filename(),this.fullname());
            members.Location=this.fullname();
        end
    end
     
end