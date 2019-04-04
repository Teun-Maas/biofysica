classdef h5file < h5group
    
    properties (Access=protected)
        filename_
        file_id_ 
    end
    
    methods
        function this = h5file(filename,options)
            if nargin < 2
                options='open';
            end
            file_id=NaN;
            if strcmpi(options,'create')
%                 if exist(filename,'file')
%                     ME=MException('h5file:h5file:fileExists','file exists: %s',filename);
%                     throw(ME);
%                 end
                fcpl = H5P.create('H5P_FILE_CREATE');
                fapl = H5P.create('H5P_FILE_ACCESS');
                file_id = H5F.create(filename, 'H5F_ACC_TRUNC', fcpl, fapl);
                
                disp(['created h5 file: ' filename]);
            elseif strcmpi(options,'open')
                if ~exist(filename,'file')
                    ME=MException('h5file:h5file:fileNotFound','file not found: %s',filename);
                    throw(ME);
                end
                file_id = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');

                disp(['opened h5 file: ' filename]);
            end
            this=this@h5group(0,'/');
            this.loc_id_=file_id;
            this.filename_=which(filename); % expand to the full pathname
        end
        
        function delete(this)
             H5F.close(this.loc_id_);
        end
        
        function fn=filename(this)
            fn=this.filename_;
        end
        
        function fid=file_id(this)
            fid=this.loc_id_;
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