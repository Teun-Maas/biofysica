% RL: This code is for RZ6 ID #3681 (4DSP) in the Vestibular Chair lab

classdef biox_vst_chr < biox_rz6_3c

    methods
        function this = biox_vst_chr
            this@biox_rz6_3c(1)          
            %RL: gemeten offsets van de DACs in [Volt]
            warning('real DACoffsets are not provided, assuming 0');
            Offsets_A(1) = 0; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) = 0; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) = 0; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) = 0; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) = 0; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) = 0; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                        
        end
    end
end
