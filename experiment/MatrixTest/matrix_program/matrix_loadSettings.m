function matrix_loadSettings(h, filename)

    s = load(filename);

    % matrix settings
    set(h.bg_listLength, 'SelectedObject', findall(h.fig_main,'tag',s.user.listLengthTag));
    set(h.bg_levelAdapt, 'SelectedObject', findall(h.fig_main,'tag',s.user.levelAdaptTag));
    set(h.cb_guess,'Value',s.user.guess);
    set(h.pupil_dilation, 'Value', s.user.pupil);
    set(h.rb_headphones, 'Value', strcmp(s.user.outputsnd, get(h.rb_headphones,'tag')));
    set(h.rb_speakers, 'Value', strcmp(s.user.outputsnd, get(h.rb_speakers,'tag')));
    try 
        set(h.edt_olsaTarget,'String',s.user.olsaTarget);
    catch ex
        set(h.edt_olsaTarget,'String','50');
    end
    try 
        set(h.edt_minStep,'String',s.user.minStep);
    catch ex
        set(h.edt_minStep,'String','0');
    end
    % calibration settings    
    set(h.edt_sensIn,'String',s.user.fullScaleOutput);
    % speech settings
    set(h.edt_speechLvl,'String',s.user.speechLevel_dB);
    set(h.edt_speechDelay,'String',s.user.speechDelay);
    try 
        set(h.cb_nomSpeechLvl,'Value',s.user.useNomSpeechLevel);
    catch ex
        set(h.cb_nomSpeechLvl,'Value',1);
    end
    % noise settings
    set(h.txt_noiseFile,'String',s.user.noiseFile);
    set(h.edt_noiseLvl,'String',s.user.noiseLevel_dB);
    set(h.bg_noiseLvlType,'SelectedObject', findall(h.fig_main,'tag',s.user.noiseLvlTypeTag));
    set(h.edt_crossDur,'String',s.user.crossDur);
    
    %channel settings
    set(h.cb_speech_L, 'Value', s.user.speechOn_L);
    set(h.cb_speech_R, 'Value', s.user.speechOn_R);
    set(h.cb_noise_L, 'Value', s.user.noiseOn_L);
    set(h.cb_noise_R, 'Value', s.user.noiseOn_R);
    set(h.cb_vocoder_L, 'Value', s.user.vocoderOn_L);
    set(h.cb_vocoder_R, 'Value', s.user.vocoderOn_R);
        
    try 
        set(h.cb_randSeekNoise,'Value',s.user.randSeekNoise);
    catch ex
        set(h.cb_randSeekNoise,'Value',0);
    end
    set(h.edt_rgb,'String',num2str(s.user.greyvalue));
    
end