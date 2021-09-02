%#########################################################################
%#                                                                       #
%#                                                                       #
%#                         matrix_prepareRZ6                            #
%#                                                                       #
%#                                                                       #
%#########################################################################

function matrix_prepareRZ6(rz6, s)

    rz6.write_wavdata([s.soundIn1, s.soundIn2], [1 2]);    

    tl = biox_rz6_tasklist;  

    tl.add_task(0.000, 'WaitForTrigger', 'Soft1'); 
    tl.add_task(0.000, 'Att',s.attenuationA, s.attenuationB);    

    % SPEECH (and Noise)
    if s.useMux
        SetMux();
    end
    
    tl.add_task(0.000, 'DAQ','Start', [1 2 3 4]);      
    tl.add_task(s.delaySND1/1000, 'SoundA', 'WAV');            
    tl.add_task(s.delaySND2/1000, 'SoundB', 'WAV');                   
    tl.add_task(0.000, 'WaitForTrigger', 'Soft2');     
    tl.add_task(0.000, 'DAQ','Stop', [1 2 3 4]); 
    
    if s.useMux
        ResetMux();
    end    
    
    tl.add_task(0.000, 'Ready');       

    rz6.write_tasklist(tl);
    delete(tl);

%#########################################################################
%#                                                                       #
%#                                                                       #
%#                          nested functions                             #
%#                                                                       #
%#                                                                       #
%######################################################################### 

    function SetMux(s)
        dev1  = s.speakers(s.indS,2);
        spk1  = s.speakers(s.indS,3);
        tl.add_task(0.000, 'MUX', dev1, 'Set', spk1);
        if s.useNoiseSpeaker
            dev2  = s.speakers(s.indN,2);
            spk2  = s.speakers(s.indN,3);
            tl.add_task(0.000, 'MUX', dev2, 'Set', spk2); 
        end  
    end 

    function ResetMux()
        tl.add_task(0.000, 'MUX', 0, 'Reset'); 
        tl.add_task(0.000, 'MUX', 1, 'Reset');
        tl.add_task(0.000, 'MUX', 2, 'Reset');           
    end

end


