% RL: This code is for RZ6 ID #4407 (4DSP) in the Speaker Arm lab


classdef biox_spk_arm < biox_rz6_3c % RL: 3c version is now identical to 4c version

    methods
        function this = biox_spk_arm
            this@biox_rz6_3c(1)
            %RL: gemeten offsets van de DACs in [Volt]
            Offsets_A(1) = -0.021;% -0.0248; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) = 0.0142;% 0.0099; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) = 0.0177;% 0.0134; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) = -0.006;% 0.0008; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) = -0.0119;%-0.0053; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) = -0.0126;%-0.0053; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                        
        end
    end
end
