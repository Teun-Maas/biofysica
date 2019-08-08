classdef rz6_tq < handle
    properties (Access=protected)
        module;
        tag_name;
        tq = [];
        nq = 0;
        ip;
    end

    properties (Constant)
        desclen = 7;

        wait_for_trigger = 0;
        start_stop_sound_a = 1;
        start_stop_sound_b = 2;
        set_mux = 3;
        set_signaling_byte = 4;
        start_stop_moving_sounds = 5;
        start_stop_daq = 6;
        set_digital_out = 7;
        output_trigger = 8;


        stop_sound = 0;
        start_tone = 1;
        start_sweep = 2;
        start_noise = 3;
        start_ripple = 4;
        start_wav = 5;
        start_ba = 6;
    end

    methods
        function this = rz6_tq(rz6_module)
            this.module = rz6_module;
            this.tag_name = 'wherethedatagoes';
        end
        
        function upload(this)
            tq_i32 = int32(this.tq(1:this.nq,:))';     
            nOS = 0;
            %TODO
            this.module.WriteTagVEX(this.tag_name, nOS, 'I32', tq_i32);
        end
        
 
        function add_task(this, delaytime, task, varargin)
           desc = this.parse_command(delaytime, task, varargin{:});
           this.append_task(desc);
        end
    end


    methods (Access=protected)
        function append_task(this, desc)
           this.nq = this.nq+1;
           this.tq(this.nq,:) = struct2array(desc);
        end

        function desc = parse_command(this, delaytime, command, varargin)
            funcName='rz6_tq/parse_command';
            validCommands = { 'WaitForTrigger', 'SoundA','SoundB','Mux','Signaling',...
               'SoundMov','Daq','SetDIO','TrigOut','Reset','Ready','HoldInp'};
            command=validatestring(command,validCommands,funcName,'Command',1);

            validateattributes(delaytime,{'numeric'},{'scalar','nonnegative'}, ...
                funcName,'DelayTime',2);

            switch lower(command)
                case 'waitfortrigger'
                    desc = this.parse_waitfortrigger(varargin{:});

                case 'sounda'
                    desc = this.parse_sounda(varargin{:});

                case 'soundb'
                    desc = this.parse_soundb(varargin{:});

                case 'mux'
                    desc = this.parse_mux(varargin{:});

                case 'signaling'
                    desc = this.parse_signaling(varargin{:});

                case 'setsignalingbyte'
                    desc = this.parse_setsignalingbyte(varargin{:});

                case 'soundmov'
                    desc = this.parse_soundmov(varargin{:});

                case 'daq'
                    desc = this.parse_daq(varargin{:});

                case 'setdio'
                    desc = this.parse_setdio(varargin{:});

                case 'trigout'
                    desc = this.parse_trigout(varargin{:});

                case 'reset'
                    desc = this.parse_reset(varargin{:});

                case 'ready'
                    desc = this.parse_ready(varargin{:});

                case 'holdinp'
                    desc = this.parse_holdinp(varargin{:});

            otherwise
                error('unknown task');
            end

            desc.DelayTime = delaytime;
            if any(isnan(struct2array(desc)))
               error('missing arguments');
            end
        end

        function desc = newdesc(this)
            desc = struct( ...
              'TaskType', 0, ...
              'SoundType', 0, ...
              'DelayTime', 0, ...
              'Par1', 0, ...
              'Par2', 0, ...
              'Par3', 0, ...
              'Par4', 0 ...
            );
        end

        function desc = parse_waitfortrigger(this,varargin)
            desc = this.newdesc();
            desc.TaskType = this.wait_for_trigger;

            p = inputParser;

            validDelayTime = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
            validExtTrig = @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative','<',8});
            %p.addRequired('DelayTime', validDelayTime);
            p.addRequired('Input', @(x) validInput(x));
            p.addOptional('ExternalTrigger', 0, @(x) validExtTrig(x));
            p.parse(varargin{:});
            %desc.DelayTime = p.Results.DelayTime;
            desc.Par1 = input2num(p.Results.Input);
            desc.Par2 = p.Results.ExternalTrigger;

            function isValid = validInput(x)
               expectedStringInputs = { 'ZBusB', 'External', 'Soft1' };
               if ischar(x)
                  isValid = any(validatestring(x,expectedStringInputs));
               else
                  isValid = validateattributes(x, {'numeric'},{'scalar','nonnegative','<',8});
               end
            end

            function n = input2num(x)
               % expectedStringInputs = { 'ZBusB', 'External', 'Soft1' };
               if ischar(x)
                  x = lower(x);
                  if x(1) == 'z'
                     n = 1;
                  elseif x(1) == 'e'
                     n = 2;
                  elseif x(1) == 's'
                     n = 4;
                  end
               else
                  n = x;
               end
            end
        end

        function desc = parse_sounda(this,varargin)
            desc = this.parse_sound(1,varargin{:});
        end

        function desc = parse_soundb(this,varargin)
            desc = this.parse_sound(2,varargin{:});
        end

        function desc = parse_sound(this,channel,varargin)
            desc = this.newdesc();
            if channel==1
               desc.TaskType = this.start_stop_sound_a;
            elseif channel==2
               desc.TaskType = this.start_stop_sound_b;
            else
               error('expected channel to be 1 or 2');
            end

            expectedSounds = { 'None','Stop','Tone','Sweep','Noise','Ripple','Wav','B=A'};

            p = inputParser;

            p.addRequired('SoundType', @(x) any(validatestring(x, expectedSounds)));
            p.addOptional('arg1', NaN);
            p.addOptional('arg2', NaN);
            p.addOptional('arg3', NaN);
            p.addOptional('arg4', NaN);
            p.parse(varargin{:});

            % Expand partially matched SoundType strings, and convert to lowercase
            SoundType=lower(validatestring(p.Results.SoundType, expectedSounds));

            switch SoundType
            case 'stop'
               desc.SoundType = 0;

            case 'tone'
               desc.SoundType = 1;
               desc.Par1 = p.Results.arg1; % f_center
               desc.Par2 = p.Results.arg2; % f_mod
               desc.Par3 = p.Results.arg3; % mod_bw
               desc.Par4 = p.Results.arg4; % att

            case 'sweep'
               desc.SoundType = 2;
               desc.Par1 = p.Results.arg1; % f_start
               desc.Par2 = p.Results.arg2; % n_octaves
               desc.Par3 = p.Results.arg3; % sweeps/1000sec
               desc.Par4 = p.Results.arg4; % att

            case 'noise'
               desc.SoundType = 3;
               desc.Par1 = p.Results.arg1; % f_lp
               desc.Par2 = p.Results.arg2; % f_hp
               desc.Par4 = p.Results.arg3; % att

            case 'ripple'
               desc.SoundType = 4;
               desc.Par1 = p.Results.arg1; % f_start
               desc.Par2 = p.Results.arg2; % mod in time
               desc.Par3 = p.Results.arg3; % mod in freq
               desc.Par4 = p.Results.arg4; % att

            case 'wav'
               desc.SoundType = 5;
               desc.Par1 = p.Results.arg1; % wav_index
               desc.Par4 = p.Results.arg2; % att

            case 'b=a'
               desc.SoundType = 6;
               desc.Par1 = p.Results.arg1; % mov_type
               desc.Par2 = p.Results.arg2; % period
               desc.Par3 = p.Results.arg3; % start_phase
               desc.Par4 = p.Results.arg4; % att

            otherwise
               error(p.Results.SoundType);
            end

        end

        function desc = parse_mux(this,varargin)
           desc = this.newdesc();
           desc.TaskType = this.set_mux;

           p = inputParser;
           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           p.addRequired('MuxByte', @(x) validByte(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.MuxByte;
        end

        function desc = parse_signaling_byte(this,varargin)
           desc = this.newdesc();
           desc.TaskType = this.set_signaling_byte;

           p = inputParser;
           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           p.addRequired('SignalingByte', @(x) validByte(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.SignalingByte;
        end

        function desc = parse_soundmov(this,varargin)
           desc = this.newdesc();
           desc.TaskType = this.start_stop_moving_sounds;

           p = inputParser;
           validBool = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<=',1});
           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           validNumSpeakers = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','>',1});
           validPeriod = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           validStartPhase = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>=',-180, '<=',180});
           p.addRequired('StartStop', @(x) validBool(x));
           p.addRequired('NumSpeakers', @(x) validNumSpeakers(x));
           p.addRequired('Period', @(x) validPeriod(x));
           p.addRequired('StartPhase', @(x) validStartPhase(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.StartStop;
           desc.Par2 = p.Results.NumSpeakers;
           desc.Par3 = p.Results.Period;
           desc.Par4 = p.Results.StartPhase;
        end

%
% !!! opsplitsen in start en stop functies
        function desc = parse_daq(this,varargin)
           desc = this.newdesc();
           desc.TaskType = this.start_stop_daq;

           p = inputParser;
           validBool = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<=',1});
           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           validAcqTime = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','positive'});
           validInputSelByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           p.addRequired('StartStop', @(x) validBool(x));
           p.addRequired('AcqTime', @(x) validAcqTime(x));
           p.addRequired('InputSelByte', @(x) validInputSelByte(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.StartStop;
           desc.Par2 = p.Results.AcqTime;
           desc.Par3 = p.Results.InputSelByte;

        end

        function desc = parse_setdio(this,varargin)
        end

        function desc = parse_trigout(this,varargin)
        end

        function desc = parse_reset(this,varargin)
        end

        function desc = parse_ready(this,varargin)
        end

        function desc = parse_holdinp(this,varargin)
        end


    end

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
