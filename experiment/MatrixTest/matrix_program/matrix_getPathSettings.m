%#########################################################################
%#                                                                       #
%#                                                                       #
%#                       matrix_getPathSettings                          #
%#                                                                       #
%#                                                                       #
%#########################################################################

% Make sure you have an OLSA directory with Settings, Results, Sounds

function path = matrix_getPathSettings
    %PATH_MATRIX = 'C:\Users\Svetlana...etc'        
    %PATH_OLSA   = 'C:\Users\Svetlana...etc'        
    path.MATRIX  = 'C:\Users\Ruurd Lof\Gitlab\matrixtest';
    path.OLSA    = 'C:\Users\Ruurd Lof\OLSA'; % RL at work 
    
    
    path.SETTINGS   = fullfile(path.OLSA, 'Settings');
    path.RESULTS    = fullfile(path.OLSA, 'Results');
    path.SOUNDS     = fullfile(path.OLSA, 'Sounds');
    path.PICTURES   = fullfile(path.OLSA, 'Pictures');
    path.LISTBASE   = fullfile(path.SOUNDS, 'Testlists_NL_20 and Testlists_NL_30');
    path.SOUNDS_NL  = fullfile(path.SOUNDS, 'Sounds_NL');
    path.CALFILE    = fullfile(path.SOUNDS_NL,  'dutchspeechnoise.wav');     
end
