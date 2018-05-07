classdef h5group < h5object
    
    properties (Access=protected)
        
    end
    
    methods
        function this = h5group(parent_group,name)
            this=this@h5object(parent_group,name);
        end
          
        function result=exists(this,link_name)
            fid = H5F.open(this.filename);
            gid = H5G.open(fid,this.fullname);
            result=H5L.exists(gid,link_name,'H5P_DEFAULT');
            H5G.close(gid);
            H5F.close(fid);
        end
        
        function group=creategroup(this,name)
            fid = H5F.open(this.filename);
            plist = 'H5P_DEFAULT';
            gid = H5G.create(fid,name,plist,plist,plist);
            H5G.close(gid);
            H5F.close(fid);
            group=h5group(this,name);
        end
        
        function group=opengroup(this,name)
            group=h5group(this,name);
        end
        
        function dataset=createdataset(this,name,varargin)
            if name(1)=='/'
                dsname=name;
            else
                dsname=strcat(this.fullname,'/',name);
            end
            h5create(this.filename,dsname,varargin{:});
            dataset=h5dataset(this,name);
        end
        
        function dataset=opendataset(this,name)
            dataset=h5dataset(this,name);
        end
        
    end
    
end