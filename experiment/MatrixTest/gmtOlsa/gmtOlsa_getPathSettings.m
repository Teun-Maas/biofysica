function path = gmtOlsa_getPathSettings
    %Root directory of all MATRIX files:        
    %PATH_MATRIX = 'C:\Users\Svetlana...etc'    
    %PATH_MATRIX = 'C:\Users\Ruurd\Documents\GitlabLocal\MATRIX_RL_2\MatlabOlsa';  %RL at home    
    path.MATRIX = 'C:\Users\Ruurd Lof\Gitlab\MATRIX_RL_2\MatlabOlsa'; % RL at work  
    
    path.SETTINGS  = fullfile(path.MATRIX,  '#Settings');
    path.RESULTS   = fullfile(path.MATRIX,  '#Results');
    path.OLSA      = fullfile(path.MATRIX,  'gmtOlsa');
    path.OLSAMAT   = fullfile(path.MATRIX,  'OlsaMaterial');        
    path.LISTBASE  = fullfile(path.OLSAMAT, 'Testlists_NL_20 and Testlists_NL_30');
    path.SOUNDS    = fullfile(path.OLSAMAT, 'Sounds_NL');
    path.CALFILE   = fullfile(path.SOUNDS,  'dutchspeechnoise.wav'); 
end
