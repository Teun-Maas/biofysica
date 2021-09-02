function matrix_ShowScreenAttr(screen, index) 

    window =screen.window;
    rect = screen.rect;    
    RGB = [screen.greyvalue screen.greyvalue screen.greyvalue];    
    Screen('FillRect', window, RGB);
    ScreenAttr = screen.attributes(index).values;    
    
    ph_MONITOR = rect(3);
    pv_MONITOR = rect(4); 
        
    baseRect  = [0 0  200  200];
    diodeRect = [0 0  350  550];
    
    if ismember('Welcome', ScreenAttr)
        Screen('TextFont', window, 'Courier');
        Screen('TextSize', window, 120);
        DrawFormattedText(window,'welkom!', 'center', 'center', 0); % 183 = dot  
    end
    
    if ismember('TopLeft', ScreenAttr)
        centeredRect = CenterRectOnPointd(baseRect,0,0); % Top left red
        Screen('FillOval',window,[255,0,0], centeredRect);
    end    
    
    if ismember('TopRight', ScreenAttr)
        centeredRect = CenterRectOnPointd(baseRect,ph_MONITOR,0); % top right
        Screen('FillOval',window,[0,0,255], centeredRect);        
    end    
    
    if ismember('BottomLeft', ScreenAttr)
        centeredRect = CenterRectOnPointd(baseRect,0,pv_MONITOR); % bottom left
        Screen('FillOval',window,[238,130,238], centeredRect);        
    end    
    
    if ismember('FixCross', ScreenAttr)
    % draw the fixation cross 
        Screen('DrawLines',window,[(ph_MONITOR/2)-100 (ph_MONITOR/2)+100 (ph_MONITOR/2) (ph_MONITOR/2); (pv_MONITOR/2) (pv_MONITOR/2) (pv_MONITOR/2)-100 (pv_MONITOR/2)+100],6, [0 0 139],[0 0]);
    end    
        
    if ismember('MatrixPic', ScreenAttr)
        paths = matrix_getPathSettings;        
        bigImSq   = [0 0 1000 1000];
        MatrixPic = fullfile(paths.PICTURES, 'Matrix.png');                     
        MatrixLoad = imread(MatrixPic); 
        MatrixImg = Screen('MakeTexture',window, MatrixLoad);
        [bigIm, ~, ~] = CenterRect(bigImSq, rect);
        Screen('DrawTexture', window, MatrixImg,[],bigIm); % draw the scene         
    end   
    
    if ismember('DiodeWhite', ScreenAttr)        
        centeredRect = CenterRectOnPointd(diodeRect,ph_MONITOR,pv_MONITOR); 
        Screen('FillRect',window,[255,255,255], centeredRect);                 
    end  
    
    if ismember('DiodeBlack', ScreenAttr)        
        centeredRect = CenterRectOnPointd(diodeRect,ph_MONITOR,pv_MONITOR); 
        Screen('FillRect',window,[0,0,0], centeredRect);                 
    end    
    
    Screen('Flip',window);        
end


