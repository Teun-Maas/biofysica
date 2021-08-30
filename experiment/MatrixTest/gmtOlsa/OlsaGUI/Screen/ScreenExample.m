function ScreenExample
  grayValue = 60;
  try            
    Screen('Preference', 'SkipSyncTests',2); % avoid timing synch
    [window, rect]=Screen('OpenWindow',2);  
    % 'Welkom'; 'TopLeft'; 'TopRight'; 'BottomLeft'; 'FixCross'; 'MatrixPic'; 'DiodeWhite'; 'DiodeBlack';   
    ScreenAttr = ['DiodeBlack', 'TopLeft', 'TopRight', 'BottomLeft', 'FixCross'];
    gmtOlsa_ShowScreen(window, rect, ScreenAttr, grayValue);    
  catch e      
    Screen('CloseAll');
    disp(getReport(e));
  end  
  WaitSecs(4);
  Screen('CloseAll');
end        