function soadata=lsl_optitrack_convert2soa(optitrack_data)
%LSL_OPTITRACK_CONVERT2SOA -  convert lsl optitrack data to structure of
%array type data

    % 16 fields
    rigid_fields= {
        'FrameID','fTimestamp', 'CameraDataReceivedTimestamp', 'TransmitTimestamp',...
        'SoftwareLatency','TransmitLatency','BodyID',...
        'Error','Valid','x','y','z','qx','qy','qz','qw'...
        };
    % 14 fields, only one timestamp
    rigid_fields_deprecated= {
        'FrameID','fTimestamp',...
        'SoftwareLatency','TransmitLatency','BodyID',...
        'Error','Valid','x','y','z','qx','qy','qz','qw'...
        };

    % 15 fields
    marker_fields= {
        'FrameID', 'fTimeStamp', 'CameraDataReceivedTimestamp', 'TransmitTimestamp',...
        'SoftwareLatency', 'TransmitLatency', 'ModelID', ...
        'MarkerID', 'Occluded', 'PCSolved', 'ModelSolved', ...
        'size', 'x', 'y', 'z'...
        };
    
    % 13 fields, only one timestamp
    marker_fields_deprecated= {
        'FrameID', 'fTimeStamp',...
        'SoftwareLatency', 'TransmitLatency', 'ModelID', ...
        'MarkerID', 'Occluded', 'PCSolved', 'ModelSolved', ...
        'size', 'x', 'y', 'z'...
        };

    if size(optitrack_data.Data,1)==numel(rigid_fields)
        fields=rigid_fields;
    elseif size(optitrack_data.Data,1)==numel(rigid_fields_deprecated)
        fields=rigid_fields_deprecated;
    elseif size(optitrack_data.Data,1)==numel(marker_fields)
        fields=marker_fields;
    elseif size(optitrack_data.Data,1)==numel(marker_fields_deprecated)
        fields=marker_fields_deprecated;
    else
        error('cannot determine OptiTrack data type');
    end

    soadata=struct();
    mfield=0;
    for fieldname=fields
        mfield=mfield+1;
        soadata.(fieldname{1})=optitrack_data.Data(mfield,:);
    end
end
