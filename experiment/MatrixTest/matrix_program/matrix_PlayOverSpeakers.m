function matrix_PlayOverSpeakers(state, rz6, settings, mixedSignal)

   msgbox('Functionality for playing over the speakers is not implemented');
   return; %skip the rest of the function
   
   silent = zeros(size(mixedSignal,1), 1);
   if (strcmp(settings.setupconf, 'rb_s0n0')) % S0N0 presentation
        signal = mixedSignal(:,1) +mixedSignal(:,2);
        signal = resample(signal,floor(rz6_sf), fsPlayback);
        if ~isrow(signal) % readTagV requires row vectors
            signal = signal';
        end                        
        nsampS=length(signal);   
        nsampN=0; %RL
        noise=0;  %RL 
        sound_duration = nsampS/fsPlayback;
        % dev and chan indices
        indS = 1;
        indN = 0;
        useNoiseSpeaker = false;
    else % spatial presentation
        signal = mixedSignal(:,1);
        signal = resample(signal,floor(rz6_sf), fsPlayback);
        useNoiseSpeaker = true;
        noise = mixedSignal(:,2);
        noise  = resample(noise,floor(rz6_sf),fsPlayback);				        
        nsampS=length(signal);
        nsampN=length(noise);
        sound_duration = nsampN/fsPlayback;
        if ~isrow(signal)
            signal = signal';
            noise = noise';
        end
        % dev and chan indices for noise depend on spatial
        % configuration
        indS = 1;
        if (strcmp(settings.setupconf, 'rb_s0nmin90')) %'S0Nmin90'
            indN = 10;
        elseif (strcmp(settings.setupconf, 'rb_s0nplus90')) %'S0Nplus90'
            indN = 26;
        elseif    (strcmp(settings.setupconf, 'rb_s0nmin70')) %'S0Nmin70'
            indN = 8;
        elseif (strcmp(settings.setupconf, 'rb_s0nplus70')) %'S0Nplus70'
            indN = 24;
        end
    end
    % Scale to +/- 10 (you will attenuate later)
    scaling = 10/max(abs(signal));
    signal = scaling*signal;
    if exist('noise','var')
        scaling = 10/max(abs(noise));
        noise = scaling*noise;
    end
		 

    attenSig=settings.maxSoundLvl-state.speechLvl(iTrial); % attenuation [dB]
    attenNoise=settings.maxSoundLvl-state.noiseLvl(iTrial);

    % prepare rz6  
    rz6settings.speakers     = settings.speakers; 

    rz6settings.useNoiseSpeaker = useNoiseSpeaker;
    rz6settings.useMux       = strcmp(settings.outputsnd, 'Speakers');

    rz6settings.delaySND1    = 0;
    rz6settings.attenuationA = attenSig;
    rz6settings.soundIn1     = signal;
    rz6settings.indS         = indS;

    rz6settings.delaySND2    = 0;
    rz6settings.attenuationB = attenNoise;
    rz6settings.soundIn2     = noise;
    rz6settings.indN         = indN;                    

    prepare_RZ6_tasklist(rz6, rz6settings);                    

    rz6.trigger('soft1'); %starts sound

    WaitSecs(sound_duration);

    rz6.trigger('soft2'); %stops 
end