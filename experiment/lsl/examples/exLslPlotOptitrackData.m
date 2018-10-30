function exLslPlotOptitrackData(optitrack_data)
    
    data=lsl_optitrack_convert2soa(optitrack_data);
    ot_ts = data.Timestamp;
    % lsl_ts = optitrack_data.Timestamps;
    lsl_ts = lsl_correct_lsl_timestamps(optitrack_data);
    
    figure(101);
    clf;
    t00=ot_ts(1);
    t10=lsl_ts(1);
    plot(lsl_ts-t10-(ot_ts-t00),'.');
    grid('on');
    
    hold('on');
    plot(data.SoftwareLatency+data.TransmitLatency,'.');
    legend('lsl\_ts-ot\_ts','SoftwareLatency+TransmitLatency')
 
    qx = data.qx;
    qy = data.qy;
    qz = data.qz;
    qw = data.qw;
    
    [az,el,rot]=quaternion2azel(qx,qy,qz,qw);
    times=lsl_ts-t10;
    figure(102);
    clf;
    plot(times, az, '.', times, el, '.', times, rot, '.');
    grid('on');
    legend('\alpha','\epsilon','\rho');
    
   % plot(times,qx,'.',times,qy,'.',times,qz,'.',times,qw,'.');
    
end

