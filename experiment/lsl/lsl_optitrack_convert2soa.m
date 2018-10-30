function soadata=lsl_optitrack_convert2soa(optitrack_data)
%LSL_OPTITRACK_CONVERT2SOA -  convert lsl optitrack data to structure of
%array type data

    rigid_fields= {
        'FrameID','fTimestamp', 'CameraDataReceivedTimestamp', 'TransmitTimestamp',...
        'SoftwareLatency','TransmitLatency','BodyID',...
        'Error','Valid','x','y','z','qx','qy','qz','qw'...
        };
    
    marker_fields= {
        'FrameID', 'fTimeStamp', 'CameraDataReceivedTimestamp', 'TransmitTimestamp',...
        'SoftwareLatency', 'TransmitLatency', 'ModelID', ...
        'MarkerID', 'Occluded', 'PCSolved', 'ModelSolved', ...
        'size', 'x', 'y', 'z'...
        };

    if size(optitrack_data.Data(1))==numel(rigid_fields)
        fields=rigid_fields;
    elseif size(optitrack_data.Data(1))==numel(marker_fields)
        fields=marker_fields;
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
