classdef h5dataset < h5object
   
    methods
        function this=h5dataset(parent_group,name)
            this=this@h5object(parent_group,name);
        end
        
        %function write(this,data,start,count,stride)
        function write(this,varargin)
            h5write(this.filename,this.fullname,varargin{:});
        end
    
        %function data=read(this,start,count,stride)
        function data=read(this,varargin)
            data=h5read(this.filename,this.fullname,varargin{:});
        end
        
    end
    
end