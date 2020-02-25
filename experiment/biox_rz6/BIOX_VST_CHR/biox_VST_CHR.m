% RL: This code is for RZ6 ID #3681 (4DSP) in the Vestibular Chair lab
% RL: This code works with 'RCX_Uni_4C_50kHz_MOV_V3.14.rcx' alone

classdef biox_VST_CHR < biox_rz6_4c_mov

    methods
        function this = biox_VST_CHR(rz6number)
            this@biox_rz6_4c_mov(rz6number)
            % create a list of [mux_id mux_index] for the speaker array.
            % sp_list must have odd nr of rows and a maximum of 21 rows
            sp_list = [2 3
                       3 3
                       2 4
                       3 5
                       2 6
                       3 7
                       2 8 
                       3 9
                       2 10
                       3 11
                       2 12
                       3 13
                       2 14
                       3 15
                       2 0
                       3 1
                       2 2];
            this.mov_spm_dac('OUT-A');     % connect middle speaker to the DAC 'OUT-A'
            this.mov_sp_list(sp_list);     % upload speaker list
            
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
