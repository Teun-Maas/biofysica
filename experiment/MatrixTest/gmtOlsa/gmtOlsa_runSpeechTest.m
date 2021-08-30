% gmtOlsa_runSpeechTest(settings)
%
% Conduct Olsa test (1 list)
%
% Change log:
%   12 Dec 2012, P.Hehrmann - created, based on K_OLSA scripts by K.-H. Dyballa
%   17 Dec 2012, PH - bug fixes, MONITOR output and plotting
%   16 Jan 2012, PH - bug fix: adjust to changes in ReadWavUnit
%   04 Mar 2013, PH - added stereo capability
%   06 Mar 2013, PH - use new OlsaMixerUnit; bug-fix: detect clipping for stereo 
%   09 Apr 2013, PH - improved backwards compatibility with Matlab (audioplayer syntax)
%   27 May 2013, PH - bug fix (crossover duration)
%   08 Aug 2013, PH - added options: initial noise seek; nominal vs. real speech RMS;
%                                    OLSA performance target + min. step size
%	05 April 2019 J. van der Heijdt - Added code for playing Olsa over rz6
%	system.
%   30 June 2021, R. Lof - Changed rz6 code to BIOX
%   25 August 2021, R. Lof - complete refactoring

%#########################################################################
%#                                                                       #
%#                                                                       #
%#                         gmtOlsa_runSpeechTest                         #
%#                                                                       #
%#                                                                       #
%#########################################################################          

function gmtOlsa_runSpeechTest(settings)  
    
    [h_Olsa_gui_matrix] = initializeMatrixGui(settings); 
    
    try

    choice = waitForUserResponse(h_Olsa_gui_matrix, 'start');     
                      
        switch choice
            case 'Start Test'
                startTest( h_Olsa_gui_matrix, settings)
            case 'Cancel'
                cancelTest(h_Olsa_gui_matrix)
        end % switch 
                
    catch e
       closeWindow(h_Olsa_gui_matrix); 
       Screen('CloseAll');
       rethrow(e);
    end %try..catch   
    
%#########################################################################
%#                                                                       #
%#                                                                       #
%#                         nested functions                              #
%#                                                                       #
%#                                                                       #
%######################################################################### 

    function startTest(h_Olsa_gui_matrix, settings)
        
        h = guihandles(h_Olsa_gui_matrix);        
        trialRecordings = [];        
        trialsCompleted = 0; 

        [mixer, rz6, state, settings, h_gui_Result, monitorInfo] = initializeTest(h, settings);        

        for iTrial = 1:settings.user.nTrials 
            
            state.iTrial = iTrial;            
            [state, Choice, trialRec] = performTrial(mixer, rz6, state, settings, h_gui_Result, h_Olsa_gui_matrix, monitorInfo);
            
            trialsCompleted = state.iTrial;
            trialRecordings{iTrial} = trialRec;
            
            switch Choice
                case 'Next'
                     state = nextTest(trialRec, state, settings);
                case 'Cancel'                     
                     cancelTest(h_Olsa_gui_matrix);
                     break;
            end %switch     
            
        end %for    

        finalizeTest(trialRecordings, settings, trialsCompleted);        
        
        function [state, Choice, trialRec] = performTrial(mixer, rz6, state, settings, h_gui_Result, h_Olsa_gui_matrix, monitorInfo)                       
            
            h = guihandles(h_Olsa_gui_matrix);
            
            printDataBeforePlay(state);

            if settings.user.MONITOR                           
                monitorInfo = displayLevelProgression(h_gui_Result,state,monitorInfo);
            end
            

            trialRec = gmtOlsa_playTrial(rz6, settings, state, mixer, h);
            Choice = waitForUserResponse(h_Olsa_gui_matrix, 'response');                                    
            trialRec.response = collectedCorrectResponses(h, settings, state);

            printDataAfterPlay(trialRec)

            updateMatrixGui(h, state, trialRec); 
            
        end %performTrial  
        
    end %startTest 

    function state = nextTest(trialRec, state, settings)
        state = gmtOlsa_updateUserState(trialRec.response.nCorrect, state, settings.parAdapt);                
    end    
        
    function cancelTest(h_Olsa_gui_matrix)
        answer = questdlg('Are you sure you want to cancel before finishing the test?',...
                   'Cancel run?', 'Yes', 'No', 'Yes');
        if strcmp(answer,'Yes') 
             closeWindow(h_Olsa_gui_matrix);
             Screen('CloseAll'); 
        end        
    end %cancelTest  
    
    function closeWindow(h_gui)
        if ishandle(h_gui)
            close(h_gui);
        end                      
    end %closeWindow           

    function choice = getMatrixGuiChoice(h)
        choice = get(h.bg_control, 'UserData');
    end  %getMatrixGuiChoice

    function updateMatrixGui(h, state, trialRec)
        hWordGroups = [h.bg_col1; h.bg_col2; h.bg_col3; h.bg_col4; h.bg_col5];        
        hOptionalInfo =  [h.txt_speechLvlLabel; h.txt_noiseLvlLabel; h.txt_correctLabel; h.txt_reversalsLabel; h.txt_totalCorrectLabel; ...
                          h.txt_speechLvl;      h.txt_noiseLvl;      h.txt_correct;      h.txt_reversals;      h.txt_totalCorrect  ];  
        set(hWordGroups, 'SelectedObject', []);
        set(hOptionalInfo, 'Visible', 'on');
        set( h.txt_sentNr,   'String', num2str(state.iTrial)); %show actual sentnumber
        set( h.txt_speechLvl,'String', sprintf('%.2f [dB]',state.speechLevel_dB)); % speech level
        set( h.txt_noiseLvl, 'String', sprintf('%.2f [dB]',state.noiseLevel_dB)); % noise level
        set( h.txt_reversals,'String', num2str(state.nReversals)); % # reversals so far
        set( h.txt_direction,'String', sprintf('%+d', state.snrDirection)); % current direction of SNR change
        set( h.txt_correct,  'String', num2str(trialRec.response.nCorrect)); %show number correct words                        
    end %updateMatrixGui

    function response = collectedCorrectResponses(h, settings, state)
        sent = cell(1,5);
        sent{1} = get(get(h.bg_col1, 'SelectedObject'),'String');
        sent{2} = get(get(h.bg_col2, 'SelectedObject'),'String');
        sent{3} = get(get(h.bg_col3, 'SelectedObject'),'String');
        sent{4} = get(get(h.bg_col4, 'SelectedObject'),'String');
        sent{5} = get(get(h.bg_col5, 'SelectedObject'),'String');
        response.file_string = { settings.sentenceList{5}{state.iTrial+settings.constants.parserRowOffset}, ...
            settings.sentenceList{6}{state.iTrial+settings.constants.parserRowOffset}, ...
            settings.sentenceList{7}{state.iTrial+settings.constants.parserRowOffset}, ...
            settings.sentenceList{8}{state.iTrial+settings.constants.parserRowOffset}, ...
            settings.sentenceList{9}{state.iTrial+settings.constants.parserRowOffset}, ...
        };
        response.nCorrect = sum(strcmp(sent, response.file_string));
        response.sent = sent;
        
        if settings.user.guess
            response = guessMissingWords(response);
        else    
            response.nCorrectGuess = 0;            
        end
        
        function response = guessMissingWords(response)
            nMissing = sum(cellfun(@isempty, response.sent)); % number of missing words
            draws = rand(1,nMissing); % uniform draws from (0,1) 
            response.nCorrectGuess = sum(draws < (1/WORDS_PER_COLUMN)); % 10% chance of success each
            response.nCorrect = response.nCorrect + response.nCorrectGuess; % add
        end %guessMissingWords
        
    end %collectCorrectResponses 

    function [mixer, rz6, state, settings, h_gui_Result, monitorInfo] = initializeTest(h, settings)
%RL!    rc = pupil_remote_control('131.174.140.72',38371);  % workaround by Gunter, check IP address in command window by running 'ping -4 dcn-pl02.local' 
        monitorInfo = [];
        if settings.user.MONITOR
            h_gui_Result = figure;
        else
            h_gui_Result = [];        
        end
        rz6 = biox_eeg_nirs;
        mixer = CreateMixerStrategy(settings);
        state = initializeState(settings); 
        set(h.pb_next, 'String', 'Next');
        set(h.pb_next, 'UserData', 'Next');
        
        if settings.user.pupil
            Screen('Preference', 'SkipSyncTests',2); % avoid timing synch
            [window, rect] =Screen('OpenWindow',1);
            settings.screen.window = window;
            settings.screen.rect = rect;
            ShowWelcomeScreen(settings.screen);
        end
                
        function state = initializeState(s)
            state.speechLevel_dB = s.user.speechLevel_dB;
            state.noiseLevel_dB = s.user.noiseLevel_dB;
            state.nReversals = 0;
            state.snrDirection = 0;
        end %initializeState
 
        function ShowWelcomeScreen(screen)
            gmtOlsa_ShowScreenAttr(screen, 1);
        end %ShowWelcomeScreen
    
    end %initializeTest

    function [h_gui] = initializeMatrixGui(settings)
       h_gui = gmtOlsa_launchMatrixGUI(settings); 
       h = guihandles(h_gui);
       hWordGroups = [h.bg_col1; h.bg_col2; h.bg_col3; h.bg_col4; h.bg_col5];        
       set(hWordGroups, 'SelectedObject', []);
       set(h.pb_next, 'String', 'Start');
       set(h.pb_next, 'UserData', 'Start Test');
       set(h.pb_replay, 'UserData', 'Replay');
       set(h.pb_replay, 'Visible', 'Off');
       set(h.pb_cancel, 'UserData', 'Cancel');
    end %initializeMatrixGui

    function finalizeTest(trialRecs, settings, trialsCompleted)
        
        recordings.trials = trialRecs;
        recordings.settings  = settings;
        recordings.trialsCompleted = trialsCompleted;
        recordings = compileRecordings(recordings); 
        saveRecordings(recordings);
        if settings.user.MONITOR
           displayResultsSummary(recordings, trialRecs);
        end   

        function recordings = compileRecordings(recordings)
            
            recordings.date               = datestr(now, 'ddmmmyyyy_HH:MM');
            recordings.speechLevelType    = 'rms (nominal)';       
            recordings.totals             = getTotalsRecording(recordings);         

            function totals = getTotalsRecording(recordings)
                SRT_nTrials = recordings.settings.constants.SRT_nTrials;
                nTrials     = recordings.trialsCompleted;
                totals.trialsSRT = max(1,(nTrials-SRT_nTrials+1)) : nTrials;                
                [~, totals.nTrialsSRT] = size(totals.trialsSRT);
                sum_SpeechLevels_dB = 0;
                sum_NoiseLevels_dB = 0;
                sum_SRT_dB = 0;
                sum_nC = 0;                
                for i = totals.trialsSRT                    
                    speechLevel_dB      = recordings.trials{i}.state.speechLevel_dB;
                    noiseLevel_dB       = recordings.trials{i}.state.noiseLevel_dB;
                    sum_SpeechLevels_dB = sum_SpeechLevels_dB + speechLevel_dB;
                    sum_NoiseLevels_dB  = sum_NoiseLevels_dB + noiseLevel_dB;   
                    sum_SRT_dB          = sum_SRT_dB + (speechLevel_dB - noiseLevel_dB);
                    sum_nC              = sum_nC + recordings.trials{i}.response.nCorrect;
                end
                totals.nCorrect = sum_nC;                
                totals.avgSpeechLevel_dB = sum_SpeechLevels_dB/totals.nTrialsSRT;
                totals.avgNoiseLevel_dB  = sum_NoiseLevels_dB/totals.nTrialsSRT; 
                totals.avgSRT_dB         = sum_SRT_dB/totals.nTrialsSRT;
            end %getTotalsRecording       

        end %compileRecordings

        function saveRecordings(recordings)

            path            = recordings.settings.path.SUBJECT;
            side_L          = recordings.settings.user.sideOn_L;
            side_R          = recordings.settings.user.sideOn_R;       
            voc_L           = recordings.settings.user.vocoderOn_L;
            voc_R           = recordings.settings.user.vocoderOn_R;
            subjectID       = recordings.settings.user.subjectId;

            fName = sprintf('%s_tr%d_side%s_voc%s_COMP_%s', subjectID, recordings.trialsCompleted,...
            num2str([side_L side_R]),num2str([voc_L, voc_R]), ...
            strrep(recordings.date,':','')); 
            fName = regexprep(fName,'[ ]','');
            fullname = fullfile(path, fName);        
            save(fullname, 'recordings');                
        end %saveRecordings
        
        function displayResultsSummary(recordings, trialRecs)
            
            h_gui_results = figure;
            try
                totals = recordings.totals;                
                description = makeDescription(recordings);
                hSub1 = subplot(7,1,1:3);
                set(hSub1, 'Visible', 'off');
                text(0, 1.1, description, 'VerticalAlignment','top','Interpreter','none');
                trialsSRT = totals.trialsSRT;
                nTrials   = recordings.settings.user.nTrials;
                subplot(7,1,4:7); hold on;
                plot(recordings.settings.user.speechLevel_dB,'ro-', 'MarkerSize', 4);
                plot(recordings.settings.user.noiseLevel_dB, 'b^-', 'MarkerSize', 4);
                line([trialsSRT(1),trialsSRT(end)], [totals.avgSpeechLevel_dB, totals.avgSpeechLevel_dB], 'Color', 'r', 'LineStyle', ':');
                line([trialsSRT(1),trialsSRT(end)], [totals.avgNoiseLevel_dB , totals.avgNoiseLevel_dB ], 'Color', 'b', 'LineStyle', ':');
                for i = trialsSRT
                    plot(trialsSRT(i), trialRecs{i}.state.speechLevel_dB,'ro');
                    plot(trialsSRT(i), trialRecs{i}.state.noiseLevel_dB, 'b^');
                end %for
                text(nTrials+1.2, totals.avgSpeechLevel_dB, sprintf('  %.1f', totals.avgSpeechLevel_dB), 'FontSize', 11);
                text(nTrials+1.2, totals.avgNoiseLevel_dB,  sprintf('  %.1f', totals.avgNoiseLevel_dB ), 'FontSize', 11);
                title('Level progression');
                xlabel('Trial #');
                ylabel('Level [dB]');
                axis tight;
                legend({'speech','noise','avg. speech', 'avg. noise'}, 'Location', 'Best');                
                drawnow;
                waitfor(h_gui_results);              
            catch e
                close(h_gui_results); 
                rethrow(e)
            end % try..catch  
            
            function desc = makeDescription(recordings)            
                ttls = recordings.totals;
                user = recordings.settings.user;
                desc = {};
                desc{end+1} = sprintf('Subject ID: %s', user.subjectId);
                desc{end+1} = sprintf('Test date: %s',  recordings.date);
                desc{end+1} = sprintf('List number: olsa%d.%02d', user.listLength, user.listNr);
                desc{end+1} = sprintf('Trials completed: %d', recordings.trialsCompleted);
                desc{end+1} = sprintf('Adaptation mode: %s', recordings.settings.parAdapt.adaptMode);
                desc{end+1} = sprintf('Words correct: %.2g%%  (%d / %d)', ttls.nCorrect/(recordings.trialsCompleted*5)*100, ttls.nCorrect, recordings.trialsCompleted*5);
                desc{end+1} = sprintf('Average speech level: %.2f dB SPL %s', ttls.avgSpeechLevel_dB, user.speechLvlTypeTag);
                desc{end+1} = sprintf('Average noise level: %.2f dB SPL %s', ttls.avgNoiseLevel_dB, user.noiseLvlTypeTag);
                desc{end+1} = sprintf('SRT: %.2f dB (%s-%s)', ttls.avgSRT_dB, user.speechLvlTypeTag, user.noiseLvlTypeTag);
            end %makeDescription            
            
        end %displayResultsSummary        

    end %finalizeTest

    function choice = waitForUserResponse(h_gui, mode)
        
        h = guihandles(h_gui);
        
        switch mode 
            case 'start'
                allWordButtons(h, 'disable');
            case 'response'
                allWordButtons(h, 'enable');       
        end %switch        
        allPushButtons(h, 'enable');
        uiwait(h_gui);    % wait for user to give response and press 'Next' or 'Cancel'          
        choice = getMatrixGuiChoice(h);
        allPushButtons(h, 'disable');
        allWordButtons(h, 'disable');        
        
        function allWordButtons(h, enable)
            hWordGroups = [h.bg_col1; h.bg_col2; h.bg_col3; h.bg_col4; h.bg_col5];  
            hWordButtons = findall(hWordGroups,'Style', 'togglebutton');
            switch enable
                case 'enable'
                    set(hWordButtons, 'enable', 'on');
                case 'disable'
                    set(hWordButtons, 'enable', 'off');
            end
        end %allWordButtons
        
        function allPushButtons(h, enable)
            hPushButtons = findall(h.bg_control, 'Style', 'pushbutton');
            switch enable
                case 'enable'
                    set(hPushButtons, 'enable', 'on');
                case 'disable'
                    set(hPushButtons, 'enable', 'off');
            end
        end   
        
    end %WaitForUserResponse

    function info = displayLevelProgression(h_gui, state, info)
        info.speechLevel_dB(state.iTrial) = state.speechLevel_dB;
        info.noiseLevel_dB(state.iTrial) = state.noiseLevel_dB; 
        figure(h_gui); hold off;
        plot(info.speechLevel_dB,'ro-');
        hold on;
        plot(info.noiseLevel_dB, 'b^-');
        title('Level progression');
        xlabel('Trial #');
        ylabel('dB');
        legend({'speech [dB]','noise [dB]'});
    end    

    function printDataBeforePlay(state)
        fprintf('Sentence #%d\n', state.iTrial)                
    end %displayMonitorDataBeforePlay

    function printDataAfterPlay(trialRec)
        sent        = cellstostr(trialRec.response.sent);
        file_string = cellstostr(trialRec.response.file_string);
        signals     = trialRec.signals;
        fprintf('     Selected: %s\n', sent)
        fprintf('      Correct: %s\n', file_string)  
        fprintf('      Digital RMS: %.2f dB (re. 1)\n', 10*log10(mean(signals.^2)));
        
        function str = cellstostr(cellarray)
            [~, Size] = size(cellarray);
            str = '';
            for i = 1:Size                
                str = sprintf('%s %s', str ,cellarray{i});
            end
        end    
                
    end %displayMonitorDataAfterPlay
           
end %gmtOlsa_runSpeechTest