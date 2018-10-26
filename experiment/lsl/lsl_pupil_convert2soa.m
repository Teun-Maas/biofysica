function soadata=lsl_optitrack_convert2soa(pupil_data)
%LSL_PUPIL_CONVERT2SOA -  convert lsl pupil-labs data to structure of
%array type data
    soadata=aos2soa(lsl_pupil_pyrepr2mat(pupil_data.Data));
end
