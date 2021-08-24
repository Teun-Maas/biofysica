classdef lsl_metadata < handle

	properties (Access=protected)
		info
	end

	methods
		function this=lsl_metadata(istream)
 			this.info=istream.inlet_info();
		end

		function delete(this)

  		end

		function as_xml(this)
			xmlstr = this.info.as_xml();
		end
	end
end
