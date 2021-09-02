function h_gui = matrix_launchMatrixGUI(s)
    h_gui = matrix_gui_matrix; %create gui
        
    figure(h_gui); %make hTest the current figure
    h = guihandles(h_gui);                    

    set(h.txt_listNr, 'String', s.user.listNr);
    set(h.txt_subjectId,'String', s.user.subjectId);
    set(h.pb_next, 'String', 'Start Test');
                     
    % disable 'MONITOR' info / labels             
    if ~s.user.MONITOR
        set( [h.txt_speechLvl, h.txt_speechLvlLabel], 'Visible', 'off'); % speech level
        set( [h.txt_noiseLvl, h.txt_noiseLvlLabel], 'Visible', 'off'); % noise level
        set( [h.txt_reversals, h.txt_reversalsLabel], 'Visible', 'off'); % # reversals so far
        set( [h.txt_direction, h.txt_directionLabel], 'Visible', 'off'); % current direction of SNR change
        set( [h.txt_correct, h.txt_correctLabel], 'Visible', 'off'); %show number correct words
        set( [h.txt_totalCorrect, h.txt_totalCorrectLabel], 'Visible', 'off'); %show number correct words
    end    
end                 