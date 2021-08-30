function gmtOlsa_PlayOverHeadphones(rz6, settings, signals)  

%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          PlayOverHeadphones                           #
%#                                                                       #
%#                                                                       #
%#########################################################################    

    
    
    gmtOlsa_prepareRZ6(rz6, getRZ6Settings(signals, settings));  

    startRZ6TaskList(rz6);

    if settings.user.pupil
        showScreen(settings,2);
%RL     startPupillometry(settings);
        WaitSecs(2);
        showScreen(settings,3);
        WaitSecs(2);
        showScreen(settings,4);        
    end %if
    
    WaitSecs(getSignalLength_sec(signals, settings));

    StopRZ6_DataAcquisition(rz6);        
    
    if settings.user.pupil        
%RL     stopPupillometry;
        WaitSecs(10);
        showScreen(settings,5); %MatrixPic
    end
    
%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          nested functions                             #
%#                                                                       #
%#                                                                       #
%#########################################################################                                  
  
    function startRZ6TaskList(rz6)
        rz6.trigger('soft1');
    end 

    function StopRZ6_DataAcquisition(rz6)
        rz6.trigger('soft2');
    end       
    
    function startPupillometry(s)
        
        rc.time_sync(0.0);
        rc.start_recording(getFileNamePupil(s));        

        function fn = getFilenamePupil(s)
            side_set    = [s.user.sideOn_L, s.user.sideOn_R];
            vocoder_set = [s.user.vocoderOn_L s.user.vocoderOn_R];
            SNR = s.user.speechLevel_dB-s.user.noiseLevel_dB;
            fname_pupil = [s.user.subjectId '_side' num2str(side_set) '_voc' num2str(vocoder_set) '_snr' num2str(SNR)];
            fn = fname_pupil(find(~isspace(fname_pupil)));
        end
        
    end %recordPupillometric

    function stopPupillometry
        rc.stop_recording;  
    end %finilizePupllometric
    
    function signal = getHeadphoneSignal(signals, s, i)
        signal = s.constants.RZ6_scaling * signals(:,i);
    end %getHeadphoneSignal

    function secs = getSignalLength_sec(signals, s)
        [size_signal_L,~] = size(signals(:,1));
        [size_signal_R,~] = size(signals(:,2));
        nsamp = max(size_signal_L, size_signal_R);    
        secs = nsamp/s.constants.fsPlayback;
    end %getSignalLength_sec 

    function rz6settings = getRZ6Settings(signals, s)
        rz6settings.speakers        = s.speakers;                     
        rz6settings.useNoiseSpeaker = false;
        rz6settings.useMux          = strcmp(s.user.outputsnd,'Speakers');                                        
        rz6settings.delaySND1       = 0;
        rz6settings.attenuationA    = 0;
        rz6settings.soundIn1        = getHeadphoneSignal(signals, s, 1);
        rz6settings.indS            = 0;                    
        rz6settings.delaySND2       = 0;
        rz6settings.attenuationB    = 0;
        rz6settings.soundIn2        = getHeadphoneSignal(signals, s, 2);
        rz6settings.indN            = 0;
    end  

    function showScreen(settings, i)        
        gmtOlsa_ShowScreenAttr(settings.screen, i);        
    end %showScreen
    
end %PlayOverHeadphones