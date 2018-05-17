function corrected_timestamps=lsl_correct_pupil_timestamps(data)
    assert(isa(data,'lsl_data'),'type mismatch: expected lsl_pupil_data or derived class.');
    assert(size(data.Data,1)>=5, 'data is probably not pupil_data');
    nsamp=size(data.Data,2);
    
    % find gaps in data, generate vector of indexes, so we can compute a vector
    % of x-axis times by  xt=indexes*deltat;
    % indexes holds the indexes of sample points that actually contain
    % data. This is useful when fitting the timestamps lateron.
    [igap,ngap]=find_gaps(data);
    indexes=generate_indexes(nsamp,igap,ngap);
    
    % generate the reference timestamps using a least squares fit followed
    % by adding an offset too keep timestamps 'causal'
    dejittered_timestamps=generate_tref_lsq(indexes,data.Timestamps);
    timecorrection=lsl_get_timecorrection(data);
    corrected_timestamps=dejittered_timestamps+timecorrection;
end


function fitted_ts=generate_tref_lsq(x,ts)
    %least squares fit of a line
    coeffs=polyfit(x,ts,1);
    fitted_ts=polyval(coeffs,x);
    % we want causality, so no early arrival of data
    % the best we can do is assuming the earliest sample is the best we
    % could achieve. This does not take into account the processing time
    % from the capturing of the camera frame until the data are available.
    [mindiff,~]=min(ts-fitted_ts);
    fitted_ts=fitted_ts+mindiff;
end

function [igap,ngap]=find_gaps(data)
    pupil_timestamps=data.Data(3,:);
    ts=pupil_timestamps-pupil_timestamps(1);
    diffs=diff(ts);
    deltat=median(diffs);
    thresh=1.5*deltat;
    igap=find(diffs>thresh);
    % gaps are in data at index igap:igap+1 (for all pos)
    ts(igap:igap+1);
    tgap=ts(igap+1)-ts(igap);
    f_ngap=tgap/deltat;
    ngap=round(f_ngap);
    % check that ngap and f_ngap are very equal, i.e. f_ngap are
    % nearly integer values.
    tolerance=0.01;
    if any(abs(ngap-f_ngap)>tolerance)
        warning('find_gaps: possibly inaccurate estimate.');
    end
end

function ind=generate_indexes(n,igap,ngap)
    increments=ones(1,n);
    increments(igap+1)=ngap;
    ind=cumsum(increments);
end


