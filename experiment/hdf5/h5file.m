classdef h5file < h5group
    
    properties (Access=protected)
        filename_
    end
    
    methods
        function this = h5file(filename,options)
            if nargin < 2
                options='open';
            end
            if strcmpi(options,'create')
                if exist(filename,'file')
                    ME=MException('h5file:h5file:fileExists','file exists: %s',filename);
                    throw(ME);
                end
                file_id = H5F.create(filename, 'H5F_ACC_EXCL', 'H5P_DEFAULT', 'H5P_DEFAULT');
                H5F.close(file_id);
                disp(['created h5 file: ' filename]);
            elseif strcmpi(options,'open')
                if ~exist(filename,'file')
                    ME=MException('h5file:h5file:fileNotFound','file not found: %s',filename);
                    throw(ME);
                end
            end
            this=this@h5group(0,'/');
            this.filename_=which(filename); % expand to the full pathname
        end
        
        function fn=filename(this)
            fn=this.filename_;
        end
        
        function group=creategroup(this,name)
            fid = H5F.open(this.filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            plist = 'H5P_DEFAULT';
            gid = H5G.create(fid,name,plist,plist,plist);
            H5G.close(gid);
            H5F.close(fid);
            group=h5group(this,name);
        end

    end
    
    methods (Access=protected)
        
        function create(this)
            file_id = H5F.create(this.filename_, 'H5F_ACC_EXCL', 'H5P_DEFAULT', 'H5P_DEFAULT');
            H5F.close(file_id);
        end
        
    end
    
end