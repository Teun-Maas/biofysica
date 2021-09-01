function mixer = CreateMixerStrategy(s)

    mixer.strategy = Strategy();    
    
    mixer.speechIn = ReadWavUnit(mixer.strategy, 'Speech_IN', 1, 1);
    mixer.noiseIn  = ReadWavUnit(mixer.strategy, 'Noise__IN', 1, 1);
    
    mixer.strategy.sp_fs = s.constants.fsPlayback; %N.B. must be set after creating speechIn and noiseIn
    
    % Left channel           
    mixer.signalOut_L = OlsaMixerUnit(mixer.strategy, 'Signal_OUT_L', getNomSpeechRms(s), s.user.noiseLvlType, ...
                   s.user.fullScaleOutput, getDelays(s), 1, s.user.crossDur, s.user.noiseSeek);     
    mixer.strategy.connectProcUnits('Speech_IN','OUTPUT_1','Signal_OUT_L','INPUT_1'); 
    mixer.strategy.connectProcUnits('Noise__IN','OUTPUT_1','Signal_OUT_L','INPUT_2');

    % Right channel           
    mixer.signalOut_R = OlsaMixerUnit(mixer.strategy, 'Signal_OUT_R', getNomSpeechRms(s), s.user.noiseLvlType, ...
                   s.user.fullScaleOutput, getDelays(s), 1, s.user.crossDur, s.user.noiseSeek);                   
    mixer.strategy.connectProcUnits('Speech_IN','OUTPUT_1','Signal_OUT_R','INPUT_1'); 
    mixer.strategy.connectProcUnits('Noise__IN','OUTPUT_1','Signal_OUT_R','INPUT_2'); 

    function level = getNomSpeechRms(s)    
        if s.user.useNomSpeechLevel
            level = std(audioread(s.path.CALFILE));
        else
            level = [];
        end
    end

    function delays = getDelays(s)           
      delays = [s.user.speechDelay 0];
    end  
               
end    