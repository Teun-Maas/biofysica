function corrected_timestamps=lsl_correct_lsl_timestamps(data)
%LSL_CORRECT_LSL_TIMESTAMPS - compute time corrected timestamps
%
% tc=lsl_correct_lsl_timestamps(data) computes time corrected timestamps
% from the information in the LSL_DATA object in data. Time correction is
% done best effort by fitting the timecorrection values to the timestamps
% in the LSL_DATA object.
%
% SEE ALSO: LSL_ESTIMATE_TIMECORRECTION, LSL_CORRECT_PUPIL_TIMESTAMPS, LSL_DATA,
% LSL_SESSION, LSL_STREAM
    timecorrection=lsl_estimate_timecorrection(data);
    corrected_timestamps=data.Timestamps+timecorrection;
end

