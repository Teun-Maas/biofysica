classdef lsl_metadata_gaze < lsl_metadata
    % lsl_metadata_gaze - simplified LSL meta data specific for pupil capture LSL
    % relay plugins >= v2.0
    % 
    % metadata = lsl_metadata_gaze(an_lsl_istream);
    % metadata has fields name, type, channel_count,
    % pupil_lsl_relay_version and a cell array of structures describing
    % each channel showing label, eye, type, unit and possibly
    % coordinate_system.
    % For a description of the fields see https://github.com/sccn/xdf/wiki/Gaze-Meta-Data
   
    properties (SetAccess=immutable)
        % these are the most important properties read from the lsl streams
        % metadata, converted to a simpler structure than the original xml
        % meta data.
       
        name
        type
        channel_count
        channel_format
        pupil_lsl_relay_version
        channel
    end
    
    methods
        function this=lsl_metadata_gaze(istream)
            this@lsl_metadata(istream);
            p=this.infostruct.info;
            this.name = p.name.Text;
            this.type = p.type.Text;
            this.channel_count = str2double(p.channel_count.Text);
            this.channel_format = p.channel_format.Text;
            this.pupil_lsl_relay_version = p.desc.pupil_lsl_relay_version.Text;
            alloc_cell = cell(1,this.channel_count);
            this.channel = struct('label',alloc_cell,'eye',alloc_cell,...
                'type',alloc_cell,'unit',alloc_cell,...
                'coordinate_system',alloc_cell);
            for count = 1:this.channel_count
                c=p.desc.channels.channel{count};
                this.channel(count).label = c.label.Text;             
                this.channel(count).eye = c.eye.Text;
                this.channel(count).type = c.type.Text;
                this.channel(count).unit = c.unit.Text;
                if isfield(c,'coordinate_system')
                    this.channel(count).coordinate_system = c.coordinate_system.Text;
                end
            end                      
        end
               
    end
end
