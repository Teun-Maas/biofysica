function soadata=lsl_optitrack_convert2soa(optitrack_data)
%LSL_OPTITRACK_CONVERT2SOA -  convert lsl optitrack data to structure of
%array type data
    fields= {
        'FrameID','Timestamp','SoftwareLatency','TransmitLatency','BodyID',...
        'Error','Valid','x','y','z','qx','qy','qz','qw'...
        };
    
    soadata=struct();
    mfield=0;
    for fieldname=fields
        mfield=mfield+1;
        soadata.(fieldname{1})=optitrack_data.Data(mfield,:);
    end
end
