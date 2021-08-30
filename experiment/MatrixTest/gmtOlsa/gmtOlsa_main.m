%#########################################################################
%#                                                                       #
%#                                                                       #
%#                            gmtOlsa_main                               #
%#                                                                       #
%#                                                                       #
%#########################################################################

% 25 August 2021, R. Lof - complete refactoring

function gmtOlsa_main
                 
    initializeMatlabEnvironment;                
    h_Olsa_gui_main = gmtOlsa_gui_main; %create and launch gui
    try         
        waitForUserResponse(h_Olsa_gui_main);   
        while(ishandle(h_Olsa_gui_main))                                  
            settings = gmtOlsa_getSettings(guihandles(h_Olsa_gui_main));
            gmtOlsa_runSpeechTest(settings); 
            waitForUserResponse(h_Olsa_gui_main);
        end %while               
    catch e
        closeWindow(h_Olsa_gui_main);
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
        paths = gmtOlsa_getPathSettings;
        addpath(genpath(paths.MATRIX));
        cd(paths.MATRIX);                
    end %initializeEnvironment     

    function waitForUserResponse(h_gui)
        
        allControls(h_gui, 'enable');
        uiwait(h_gui);                        
        allControls(h_gui, 'disable');        
        
        function allControls(h_gui, enable)
            switch enable
                case 'enable'  
                    set(findall(h_gui, 'type', 'uicontrol'), 'Enable', 'on');    
                case 'disable' 
                    set(findall(h_gui, 'type', 'uicontrol'), 'Enable', 'off');  
            end %switch
        end %allControls
               
    end  %WaitForUserResponse          
     
end %gmtOlsa_main



