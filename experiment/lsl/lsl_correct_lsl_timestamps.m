function corrected_timestamps=lsl_correct_lsl_timestamps(data)
% LSL_CORRECT_TIMESTAMPS    
    timecorrection=lsl_get_timecorrection(data);
    corrected_timestamps=data.Timestamps+timecorrection;
end

 