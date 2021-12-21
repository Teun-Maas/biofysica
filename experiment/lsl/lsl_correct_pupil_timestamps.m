function corrected_timestamps=lsl_correct_pupil_timestamps(data)
    %    LSL_CORRECT_PUPIL_TIMESTAMPS - compute time corrected timestamps for
    %    Pupil Labs data acquired with the pupil_lsl plugin in PupilCapture
    %
    %    LSL_CORRECT_PUPIL_TIMESTAMPS uses LSL_CORRECT_LSL_TIMESTAMPS for
    %    correction time delay and drift of the source data.
    %The PupilLabs data however has two other problems to be taken care of:
    %1) PupilLabs data sometimes contains gaps when processing of the raw
    %camera images takes too long.
    %2) LSL timestamps generated by the pupil_lsl plugin in PupilCapture
    %contain a lot of jitter. To eliminate the jitter a least-squares fit is
    %done to find a monotonically rising time stamp series with equidistant
    %values, except for the positions where gaps occur. It is taken care for
    %that the fitted time stamps are never earlier than expected in a causal
    %system.
    %
    % tc=lsl_correct_pupil_timestamps(data) returns the best effort correction
    % for the available timing information in the LSL_DATA in DATA.
    %
    % SEE ALSO: LSL_ESTIMATE_TIMECORRECTION, LSL_CORRECT_LSL_TIMESTAMPS, LSL_DATA,
    % LSL_SESSION, LSL_STREAM
    
    %20211215 GW Changed identification of data structure, either ancient
    %            5-row matrix, the successor that is an array of structures
    %            or the pupil lsl plugin v2.1 22-row matrix in XDF format.
    
    %20211216 GW Added warning for negative jumps in timestamps
    %            This is probably a perfomance issue with Pupil Capture 3.x
    %            on older computers.
    
    %20211221 GW Improved robustness of estimation algorithm
    %            Slightly better handling of negative jumps in timestamps
    
    assert(isa(data,'lsl_data'),'type mismatch: expected lsl_pupil_data or derived class.');
    if isstruct(data.Data)
        assert(isfield(data.Data,'topic'));
    elseif (iscell(data.Data) && (data.Data{1}(1)=='{'))
        error('data may be pupil_data in Python string representation, convert using lsl_pupil_convert2soa first');
    else
        assert(size(data.Data,1)>=5, 'data is probably not pupil_data');
    end
    if isstruct(data.Data)
        nsamp=size(data.Data.timestamp,1);
    else
        nsamp=size(data.Data,2);
    end
    % find gaps in data, generate vector of indexes, so we can compute a vector
    % of x-axis times by  xt=indexes*deltat;
    % indexes holds the indexes of sample points that actually contain
    % data. This is useful when fitting the timestamps lateron.
    [igap,ngap]=find_gaps(data);
    indexes=generate_indexes(nsamp,igap,ngap);
    
    % generate the reference timestamps using a least squares fit followed
    % by adding an offset too keep timestamps 'causal'
    dejittered_timestamps=estimate_timestamps_lsq(indexes,data.Timestamps);

    timecorrection=lsl_estimate_timecorrection(data);
    corrected_timestamps=dejittered_timestamps+timecorrection;
    
%     sf=gcf;
%     figure(101);
%     plot(indexes,corrected_timestamps,'.'); grid on;
%     figure(sf);
end


function fitted_ts=estimate_timestamps_lsq(x,ts)
    dts=diff(ts);
    negative_dts = (dts<0);
    if ~any(negative_dts)
        % we have monotonically increasing timestamps, all well
        %least squares fit of a line
        coeffs=polyfit(x,ts,1);
        fitted_ts=polyval(coeffs,x);
        % we want causality, i.e. no early arrival of data.
        % The best we can do is assuming the earliest sample is the best we
        % could achieve. This does not take into account the processing time
        % from the capturing of the camera frame until the data are available.
        [mindiff,~]=min(ts-fitted_ts);
        fitted_ts=fitted_ts+mindiff;
    else
        % we have a problem here, let's do our best
        warning('lsl_correct_pupil_timestamps:estimate_timestamps_lsq: %d negative delta_t''s found in timestamp series. Probably inaccurate estimate.',sum(dts<0));
        % we don't have monotonically increasing timestamps, that is not okay
        % we can do the fit, as above, but the causality trick will leave us with a nonrealistic estimate

        %least squares fit of a line
        coeffs=polyfit(x,ts,1);
        fitted_ts=polyval(coeffs,x);
        % we want causality, i.e. no early arrival of data.
        % The best we can do is assuming the earliest sample is the best we
        % could achieve. This does not take into account the processing time
        % from the capturing of the camera frame until the data are available.
        %   [mindiff,~]=min(ts-fitted_ts); % This will give us a way of
        %   estimate because of the negative time jump. Let's leave it at 0
        mindiff=0;
        fitted_ts=fitted_ts+mindiff;

    end
    
end

function [igap,ngap]=find_gaps(data)
    if isstruct(data.Data)                      % v1 plugin, python representation format converted to SOA
        pupil_timestamps=data.Data.timestamp;   
    elseif size(data.Data,1)==5                 % first generation LSL relay, 5xN matrix
        pupil_timestamps=data.Data(3,:);
    elseif size(data.Data,1)==22                % v2.1 LSL relay plugin
        pupil_timestamps=data.Timestamps;
    else
        error('unknown pupil data format');
    end
    ts=pupil_timestamps-pupil_timestamps(1);
    diffs=diff(ts);
    if any(diffs<0)
        warning('lsl_correct_pupil_timestamps:find_gaps: %d negative delta_t''s found in timestamp series.',sum(diffs<0));
    end
    
    deltat=median(diffs);
    thresh=1.5*deltat;
    igap1=find(diffs>thresh);

    %ommit gaps that are followed by a very short deltat, i.e. sampling
    %caught up again
    if igap1(end)+2 > length(ts)
       igap1 = igap1(1:end-1);  % prevent indexing error below
    end
    tgap2=ts(igap1+2)-ts(igap1);
    thresh2=2.5*deltat;
    igap=find(tgap2 > thresh2);

    % gaps are in data at index igap:igap+1 (for all pos)
    tgap=ts(igap+1)-ts(igap);

    f_ngap=tgap/deltat;
    ngap=round(f_ngap);
    % check that ngap and f_ngap are very equal, i.e. f_ngap are
    % nearly integer values.
    tolerance=0.01;
    if any(abs(ngap-f_ngap)>tolerance)
        warning('lsl_correct_pupil_timestamps:find_gaps: possibly inaccurate estimate.');
    end
end

function ind=generate_indexes(n,igap,ngap)
    increments=ones(1,n);
    increments(igap+1)=ngap;
    ind=cumsum(increments);
end


