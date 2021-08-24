classdef lsl_gaze_metadata < lsl_metadata

	methods
		function this=lsl_gaze_metadata(istream)
			this@lsl_metadata(istream);
			fprintf('The channel labels are as follows:\n');
			ch = this.info.desc().child('channels').child('channel');
			for k = 1:this.info.channel_count()
    			fprintf(['  ' ch.child_value('label') '\n']);
    			ch = ch.next_sibling();
		end

	end


end
