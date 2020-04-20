% RL: This code is for RZ6 ID #3960 (3DSP) in the Patient lab

classdef biox_pat_lab < biox_rz6_3c

    methods
        function this = biox_PAT_LAB
            this@biox_rz6_3c(1)
            %RL: gemeten offsets van de DACs in [Volt]
            Offsets_A(1) = -0.0001; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) = -0.0040; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) = -0.0044; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) = -0.0010; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) =  0.0009; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) =  0.0105; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                        
        end
    end
end
