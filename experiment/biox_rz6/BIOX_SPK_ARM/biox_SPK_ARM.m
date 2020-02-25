% RL: This code is for RZ6 ID #4407 (4DSP) in the Speaker Arm lab
% RL: This code works with 'RCX_Uni_3C_50kHz_V3.14.rcx' alone

classdef biox_SPK_ARM < biox_rz6_4c

    methods
        function this = biox_SPK_ARM(rz6number)
            this@biox_rz6_4c(rz6number)
            %RL: gemeten offsets van de DACs in [Volt]
            Offsets_A(1) = -0.0248; %RL: offset DAC-A bij 0  dB AttA
            Offsets_A(2) =  0.0099; %RL: offset DAC-A bij 20 dB AttA
            Offsets_A(3) =  0.0134; %RL: offset DAC-A bij 40 dB AttA
            Offsets_B(1) =  0.0008; %RL: offset DAC-B bij 0  dB AttB
            Offsets_B(2) = -0.0053; %RL: offset DAC-B bij 20 dB AttB
            Offsets_B(3) = -0.0053; %RL: offset DAC-B bij 40 dB AttB
            this.write_DACoffsets(Offsets_A, Offsets_B);                        
        end
    end
end