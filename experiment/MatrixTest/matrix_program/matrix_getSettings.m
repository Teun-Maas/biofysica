function s = matrix_getSettings(guiHandles)

%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          gmtOlsa_getSettings                          #
%#                                                                       #
%#                                                                       #
%#########################################################################

    s.path         = matrix_getPathSettings;
    s.constants    = getConstants;
    s.user         = getUserSettings(guiHandles, s);
    s.parAdapt     = getParAdapt(s);
    s.sentenceList = getSentenceList(s);
    s.speakers     = getSpeakerSettings(s);
    s.path.SUBJECT = makeSubjectPath(s);
    s.screen       = getScreen(s);
    
%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          nested functions                             #
%#                                                                       #
%#                                                                       #
%#########################################################################        
       
    function constants = getConstants
        constants.fsPlayback            = 48830;    
        % Recalibrate noise (GWN) to speech levels [-7.5 dB].
        % This was done after finding during output calibration that
        % speech was still 7.5 dB softer than noise.           
        constants.noiseCorrection_dB    = -7.5;   
        constants.RZ6_scaling           = 10; %RZ6 max output is 10V  
        constants.parserRowOffset       = 4;  %skip 4 rows in sentence file
        constants.SRT_nTrials           = 10; %max nr of trials for level adaptation
        constants.maxSoundLvl           = 75; %RL was nowhere specified is used for speakers only;  
        constants.nVocoderChannels      = 13; %RL maximum value is 13 
    end %getConstants

    function user = getUserSettings(h, s)                
                
        user.subjectId          = get(h.edt_subjectId, 'String');

        user.MONITOR            = get(h.cb_monitor,'Value');    

        %list settings
        user.listNr             = get(h.pop_listNr, 'Value');        
        user.listLengthTag      = get(get(h.bg_listLength, 'SelectedObject'), 'tag');        
        
        user.levelAdaptTag      = get(get(h.bg_levelAdapt, 'SelectedObject'), 'tag');
        user.guess              = get(h.cb_guess,'Value');
        user.minStep            = str2double(get(h.edt_minStep,'String'));
        user.olsaTarget         = str2double(get(h.edt_olsaTarget,'String'));
        
        % calibration settings
        user.fullScaleOutput    = str2double(get(h.edt_sensIn,'String'));
        
        % speech settings
        user.speechLevel_dB     = str2double(get(h.edt_speechLvl,'String'));  
        user.speechDelay        = str2double(get(h.edt_speechDelay,'String'));
        user.useNomSpeechLevel  = get(h.cb_nomSpeechLvl,'Value');   
        user.speechLvlTypeTag   = 'rms (nominal)'; %no selection field implemented

        % noise settings                 
        filename                = get(h.txt_noiseFile,'String');
        [~, nFname, nFext]      = fileparts(filename);         
        user.noiseFile          = strcat(s.path.SOUNDS_NL,'\', nFname, nFext);
        user.noiseLevel_dB      = str2double(get(h.edt_noiseLvl,'String'));    
        user.noiseLvlTypeTag    = get(get(h.bg_noiseLvlType,'SelectedObject'), 'tag');
        user.crossDur           = str2double(get(h.edt_crossDur,'String'));
        user.randSeekNoise      = get(h.cb_randSeekNoise, 'Value');
        user.setupconf          = get(get(h.setupconf,'SelectedObject'), 'tag');
        user.outputsnd          = get(get(h.outputsound,'SelectedObject'), 'tag');

        % HEADPHONES Side
        user.speechOn_L         = get(h.cb_speech_L,'Value');
        user.speechOn_R         = get(h.cb_speech_R,'Value');
        user.noiseOn_L          = get(h.cb_noise_L,'Value');
        user.noiseOn_R          = get(h.cb_noise_R,'Value');
        user.vocoderOn_L        = get(h.cb_vocoder_L,'Value');
        user.vocoderOn_R        = get(h.cb_vocoder_R,'Value');

        user.pupil              = get(h.pupil_dilation,'Value');

        user.greyvalue          = str2double(get(h.edt_rgb,'String'));        

        switch user.listLengthTag
            case 'rb_list10'
                user.nTrials = 10;
            case 'rb_list20'
                user.nTrials = 20;
            case 'rb_list30'
                user.nTrials = 30;
            case 'rb_list50'
                user.nTrials = 50;
            otherwise
                error('Selected list length not supported.');
        end %switch
        
        user.listLength      = user.nTrials;

        if user.randSeekNoise
            user.noiseSeek = -1; % random
        else
            user.noiseSeek = 0;  % none
        end %if else

            % noise level type
        switch user.noiseLvlTypeTag
            case {'rb_noiseLvlRms'}
                user.noiseLvlType = 'rms';
            case {'rb_noiseLvlPeak'}
                user.noiseLvlType = 'peak';
            otherwise
                error('Selected noise level type not supported.');
        end %switch

        switch user.levelAdaptTag
            case{'rb_adaptOff'}
                user.levelAdaptMode = 'off';
            case{'rb_adaptSpeech'}
                user.levelAdaptMode = 'speech';
            case{'rb_adaptNoise'}
                user.levelAdaptMode = 'noise';
            otherwise
                error('Invalid adaptation mode: %s', user.levelAdaptTag);
        end %switch

        user.speechDelay        = str2double(get(h.edt_speechDelay,'String'));
        user.speechDelayRange   = [user.speechDelay user.speechDelay];  
 

        % for sides and vocoders there are only 6 valid options out of 
        % option 1: side = [1 0] voc = [0 0]    
        % option 2: side = [0 1] voc = [0 0]    
        % option 3: side = [1 0] voc = [1 0]    
        % option 4: side = [0 1] voc = [0 1]
        % option 5: side = [1 1] voc = [1 0]    
        % option 5: side = [1 1] voc = [0 1]    
        % option 6: side = [1 1] voc = [1 1]

        assert(max([user.speechOn_L, user.speechOn_R]) > 0, 'At least one side should be selected');
        assert(min([user.speechOn_L, user.speechOn_R] - [user.vocoderOn_L, user.vocoderOn_R]) >= 0, 'The combination of sides and vocoders is not valid'); 

    end %addUserSettings    
        
    function parAdapt = getParAdapt(settings)
        parAdapt.target     = settings.user.olsaTarget/100;
        parAdapt.slope      = 0.15;
        parAdapt.fi_scale   = 1.5;
        parAdapt.fi_base    = 1.4;
        parAdapt.minStep    = settings.user.minStep;
        parAdapt.adaptMode  = settings.user.levelAdaptMode;   
    end  %addParAdaptSettings      
    
    function list = getSentenceList(settings)
        listFile = fullfile(sprintf('%s', settings.path.LISTBASE),...
                            sprintf('nlmatrix%d.%02d', settings.user.nTrials, settings.user.listNr) );
        fid = fopen( listFile, 'r' );
        list = textscan(fid, '%s %s %s %s %s %s %s %s %s', 'Delimiter', ':/ ');
        fclose(fid);   
    end %getSentenceList

    function speakers = getSpeakerSettings(settings)
        if (strcmp(settings.user.outputsnd, 'rb_speakers'))
           cfg = spherelookupMinor();    
           speakers = cfg.lookup;
        else
           speakers = [];
        end %if else       
    end % getSpeakerSettings

    function path = makeSubjectPath(settings)
        path = fullfile(settings.path.RESULTS, settings.user.subjectId);
        if ~exist(path, 'dir')
           mkdir(path);
        end %if   
    end  %makeSubjectPath

    function screen = getScreen(settings)
        screen.attributes(1).values = ['Welcome'];
        screen.attributes(2).values = ['FixCross'];
        screen.attributes(3).values = ['DiodeBlack'];
        screen.attributes(4).values = ['DiodeWhite', 'TopLeft', 'TopRight', 'BottomLeft', 'FixCross'];
        screen.attributes(5).values = ['MatrixPic'];
        screen.greyvalue            = settings.user.greyvalue;
        settings.screen.window      = [];
        settings.screen.rect        = [];
    end    
    
end %gmtOlsa_getSettings
