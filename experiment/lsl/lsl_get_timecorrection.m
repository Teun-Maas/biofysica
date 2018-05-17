function timecorrection=lsl_get_timecorrection(data)
    assert(isa(data,'lsl_data'),'type mismatch: expected lsl_data or derived class.');
    
    x=data.TCindex;
    times=data.Timestamps(x);
    
    coeffs=polyfit(times,data.TimeCorrection,1);
    %%% CHECK HERE IF coeffs(1) is sufficiently small!!!
    %%% This happens when the number of sampled timecorrection values is too small,
    %%% and computed clock drift may be unrealistic due to a bad fit.
    %%% ->fallback to some other algorithm (use the average?)
    %%% ->or take care of having more timecorrection values.
    if coeffs(1) > 1e-4
        warning('timecorrection clock drift may be unrealistic');
    end
    timecorrection=polyval(coeffs,data.Timestamps);
end