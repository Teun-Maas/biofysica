function ScreenTest
  grayValue = 60;
  try            
    Screen('Preference', 'SkipSyncTests',2); % avoid timing synch
    [window, rect]=Screen('OpenWindow',2);  
    % 'Welkom'; 'TopLeft'; 'TopRight'; 'BottomLeft'; 'FixCross'; 'MatrixPic'; 'DiodeWhite'; 'DiodeBlack';   
    ScreenAttr = ['DiodeWhite', 'TopLeft', 'FixCross'];
    gmtOlsa_ShowScreen(window, rect, ScreenAttr, grayValue);    
  catch e      
    Screen('CloseAll');
    disp(getReport(e));
  end  
  WaitSecs(2);
  Screen('CloseAll');
end        