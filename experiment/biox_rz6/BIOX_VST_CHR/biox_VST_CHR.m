% RL: This code is for RZ6 ID #3681 (4DSP) in the Vestibular Chair lab

classdef biox_VST_CHR < biox_rz6_4c_mov

    methods
        function this = biox_VST_CHR
            this@biox_rz6_4c_mov(1)
            % create a list of [mux_id mux_index] for the speaker array.
            % sp_list must have odd nr of rows and a maximum of 21 rows
            sp_list = [ 0 0
                        1 0
                        0 1
                        1 1
                        0 2
                        1 2
                        0 3
                        1 3
                        0 4
                        1 4
                        0 5
                        1 5
                        0 6 
                        1 6
                        0 7
                        1 7
                        0 8
                        1 8
                        0 9
                        1 9
                        0 10];
                    
            this.mov_spm_dac('OUT-A');% connect middle speaker to the DAC 'OUT-A'
            this.mov_sp_list(sp_list);% upload speaker list
            
            %RL: gemeten offsets van de DACs in [Volt]
            Offsets_A(1) = ; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) = ; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) = ; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) = ; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) = ; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) = ; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                        
        end
    end
end
