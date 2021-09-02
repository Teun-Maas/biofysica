% open matrix_getPathSettings first to check the paths

%#########################################################################
%#                                                                       #
%#                                                                       #
%#                            matrix_main                                #
%#                                                                       #
%#                                                                       #
%#########################################################################

% 25 August 2021, R. Lof

function matrix_main
                 
    initializeMatlabEnvironment;                
    h_gui_main = matrix_gui_main; %create and launch main gui
    try         
        waitForUserResponse(h_gui_main);   
        while(ishandle(h_gui_main))                                  
            settings = matrix_getSettings(guihandles(h_gui_main));
            matrix_runSpeechTest(settings); 
            waitForUserResponse(h_gui_main);
        end %while               
    catch e
        closeWindow(h_gui_main);
        disp(getReport(e));        
    end %try...catch
    
%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          nested functions                             #
%#                                                                       #
%#                                                                       #
%#########################################################################    
    

    function closeWindow(h)            
        if ishandle(h)
            close(h);
        end                      
    end %closeWindow

    function initializeMatlabEnvironment
        clc;
        paths = matrix_getPathSettings;
        addpath(genpath(paths.MATRIX));
        cd(paths.MATRIX);                
    end %initializeEnvironment     

    function waitForUserResponse(h_gui)
        
        allControls(h_gui, 'enable');
        uiwait(h_gui);                        
        allControls(h_gui, 'disable'); 
        
        %##################nested functions#########################
        
        function allControls(h_gui, enable)
            switch enable
                case 'enable'  
                    set(findall(h_gui, 'type', 'uicontrol'), 'Enable', 'on');    
                case 'disable' 
                    set(findall(h_gui, 'type', 'uicontrol'), 'Enable', 'off');  
            end %switch
        end %allControls
               
    end  %WaitForUserResponse          
     
end %matrix_main



