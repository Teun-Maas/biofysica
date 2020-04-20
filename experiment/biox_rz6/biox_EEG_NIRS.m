% RL: This code is for RZ6 ID #2170 (3DSP) in the EEG/NIRS lab

classdef biox_eeg_nirs < biox_rz6_3c

    methods
        function this = biox_eeg_nirs             
            this@biox_rz6_3c(1)
            %RL: gemeten offsets van de DACs in [Volt]
            Offsets_A(1) = -0.0119; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) =  0.0119; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) =  0.0143; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) = -0.0007; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) =  0.0019; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) =  0.0021; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                  
        end
    end
end
