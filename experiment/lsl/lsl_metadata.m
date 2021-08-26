classdef lsl_metadata < handle
% lsl_metadata - Generic parser for LSL stream meta data 
% https://github.com/labstreaminglayer/liblsl-Matlab/tree/master/examples
	properties (Access=protected)
		infoxml
        infostruct
	end

	methods
		function this=lsl_metadata(istream)
 			streaminfo=istream.info();
            this.infostruct = lsl_streaminfo2struct(streaminfo);
            this.infoxml = streaminfo.as_xml();
		end

		function xmlstr = as_xml(this)
			xmlstr = this.infoxml;
        end
        
        function infostruct = as_struct(this)
            infostruct = this.infostruct;
        end
	end
end
