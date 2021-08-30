%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          gmtOlsa_playTrial                            #
%#                                                                       #
%#                                                                       #
%#########################################################################

%   25 August 2021, R. Lof - complete refactoring


function trialRec = gmtOlsa_playTrial(rz6, settings, state, mixer, h)
   
    signals = calcMixerSignals(mixer, state, settings);
    signals = calcVocoderSignals(signals, settings);
    showSignalClipping(h, signals);           
                     
    if (strcmp(settings.user.outputsnd, 'rb_headphones'))
       gmtOlsa_PlayOverHeadphones(rz6, settings, signals);       
       %RL for now SPEAKERS not availlable       
      elseif (strcmp(settings.user.outputsnd,'Speakers'))
          gmtOlsa_PlayOverSpeakers(rz6, settings, signals);
    end
       
    trialRec = collectTrialRecordings(rz6, state, signals);    
    
    
%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          nested functions                             #
%#                                                                       #
%#                                                                       #
%#########################################################################    

    function mixerSignals = calcMixerSignals(mixer, state, settings)  
        
        speechLevel_L_dB = settings.user.sideOn_L*state.speechLevel_dB;
        speechLevel_R_dB = settings.user.sideOn_R*state.speechLevel_dB;
        noiseLevel_L_dB  = settings.user.sideOn_L*state.noiseLevel_dB + settings.constants.noiseCorrection_dB;
        noiseLevel_R_dB  = state.noiseLevel_dB + settings.constants.noiseCorrection_dB;                
        
        mixer.signalOut_L.lvlDb = [speechLevel_L_dB, noiseLevel_L_dB];            
        mixer.signalOut_R.lvlDb = [speechLevel_R_dB, noiseLevel_R_dB];
        
        mixer.strategy.resetDataUnits();   
        
        sentenceFile = getSentenceFile(state, settings);        
        
        mixer.speechIn.setData('INPUT_1', sentenceFile);
        mixer.noiseIn.setData( 'INPUT_1', settings.user.noiseFile);    
                
        mixer.strategy.run();

        signal_L = mixer.signalOut_L.getData('OUTPUT_1');
        signal_R = mixer.signalOut_R.getData('OUTPUT_1');

        mixerSignals = [signal_L signal_R];  
        
        function file = getSentenceFile(state, settings)        
            row  = state.iTrial + settings.constants.parserRowOffset;
            file = fullfile(settings.path.SOUNDS, strtrim(settings.sentenceList{2}{row}));
        end %getSentenceFile
        
    end %calcMixerSignals

    function soundOut = rescale_snd(SoundIn, dB)
       scalingF    = 10^(dB/20); 
       soundOut = scalingF * SoundIn;
    end %rescale_snd

    function signalsOut = calcVocoderSignals(signalsIn, settings)
        
        signalsOut = signalsIn; % default
        
        if settings.user.vocoderOn_L  
            signalsOut(1) = calcVocoder(signalsIn(1)); % 1 = left
        else
            signalsOut(1) = signalsIn(1);
        end %if else
        
        if settings.user.vocoderOn_R  
            signalsOut(2) = calcVocoder(signalsIn(2)); % 2 = right            
        else
            signalsOut(2) = signalsIn(2);
        end %if else        
        
        function voc = calcVocoder(signal)        
            nchannels = 6;
            voc = pa_vocoderCISim(signal,'nChans',nchannels);
            voc = voc(1:size(speech,1)); % if vocoder is always longer than speech...                
            speech_dB = 20*log10(rms(signal));
            vocoder_dB = 20*log10(rms(voc));      
            voc = rescale_snd(voc, speech_dB - vocoder_dB);
        end %calcVocoder
                
    end %calcVocoderSignals

    function showSignalClipping(h, signals)               
        if any(abs(signals(:)) >= 1)
            set(h.txt_clipOut, 'Visible', 'on');    
        else
            set(h.txt_clipOut, 'Visible', 'off');
        end   
    end %showSignalClipping  

    function trialRec = collectTrialRecordings(rz6, state, signals)                        
        trialRec.acqData = rz6.read_acqdata([1 2 3 4]);                 
        trialRec.state  = state;
        trialRec.signals = signals;
    end %collectRecordings
        
end