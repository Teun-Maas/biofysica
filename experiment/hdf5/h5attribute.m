classdef h5attribute < h5location
      
    methods
        function this = h5attribute(parent, name, value)
            % location must be a group or a dataset or named datatype
            % Attributes are small, and stored in the object header of the
            % object it describes and thus the location points to this
            % object, and name identifies the attribute.
            this=this@h5location(parent, name);        
            if nargin >= 3
                this.write(value);
            end
        end
        
        function attval=read(this)
            attval=h5readatt(this.filename,this.location,this.name);
        end
        
        function write(this,attval)
            h5writeatt(this.filename,this.location,this.name,attval);
        end
        
        function p=location(this)
            p=this.parent_location.fullname();
        end
        
        
        function inf=info(this)
            inf=struct(...
                'Filename',this.filename,...
                'Name',this.name,...         
                'Type','Attribute',...
                'Location',this.location,...
                'Value',{this.read} );
        end
    end
    
end