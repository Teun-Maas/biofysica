function timecorrection=lsl_estimate_timecorrection(data)
% LSL_GET_TIMECORRECTION - compute time correction values from data
%
% tc=lsl_get_timecorrection(data) computes time corrections values from the
% LSL_DATA in data. This is done by least squares fitting the
% time correction information to the time stamps in data.
%
% tc=lsl_get_timecorrection(data) get time corrections from data. data must
% be an LSL_DATA object or derived class.
%
%SEE ALSO:  LSL_CORRECT_LSL_TIMESTAMPS,
%LSL_CORRECT_PUPIL_TIMESTAMPS, LSL_DATA, LSL_SESSION, LSL_STREAM

    assert(isa(data,'lsl_data'),'type mismatch: expected lsl_data or derived class.');
    
    x=data.TCindex;
    times=data.Timestamps(x)-data.Timestamps(1);
    coeffs=polyfit(times,data.TimeCorrection,1);
    %%% CHECK HERE IF coeffs(1) is sufficiently small!!!
    %%% This happens when the number of sampled timecorrection values is too small,
    %%% and computed clock drift may be unrealistic due to a bad fit.
    %%% ->fallback to some other algorithm (use the average?)
    %%% ->or take care of having more timecorrection values.
    if coeffs(1) > 1e-4
        warning('timecorrection clock drift may be unrealistic');
    end
    timecorrection=polyval(coeffs,data.Timestamps-data.Timestamps(1));
end