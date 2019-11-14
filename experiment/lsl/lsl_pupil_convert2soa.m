function soadata=lsl_pupil_convert2soa(pupil_data)
%LSL_PUPIL_CONVERT2SOA -  convert lsl pupil-labs data to structure of
%array type data
    soadata=aos2soa(lsl_pupil_pyrepr2mat(pupil_data.Data));
end

%function soadata_do_not_use=lsl_pupil_convert2soa(pupil_data)
%%LSL_PUPIL_CONVERT2SOA -  convert lsl pupil-labs data to structure of
%%array type data
%    if nargout > 0
%       error('sorry, calling syntax changed');
%    end
%    %%FIXME pupil_data.pythonData=pupil_data.Data;
%    pupil_data.Data = aos2soa(lsl_pupil_pyrepr2mat(pupil_data.Data));
%end
