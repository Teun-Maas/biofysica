classdef h5location < handle
    
    properties (Access = protected)
        parent_location h5location
        name_
    end
    
    methods
        function this=h5location(parent_location,name)
            if isa(parent_location,'h5location')
                this.parent_location=parent_location;
            elseif parent_location==0
                % called from h5file constructor
                %FIXME - add a dummy location or something maybe?
            else
                error('parent_location must be a location object');
            end
            this.name_=name;
        end
        
        function fn=filename(this)
            fn=this.parent_location.filename();
        end
        
        function fid=file_id(this)
            fid=this.parent_location.file_id();
        end
        
        function n=name(this)
            n=this.name_;
        end
        
        function p=fullname(this)
            if this.name_(1)=='/'
                p=this.name_;
            else
                pfn=this.parent_location.fullname();
                if pfn=='/'
                    pfn='';
                end
                p=strcat(pfn,'/',this.name_);
            end
        end
        
        function e=exists(this,name)
            
        end
    end
end