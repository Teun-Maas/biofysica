classdef rz6_unircx_tasklist < handle
    properties (Access=protected)
        module;
        tag_name;
        tq = [];
        nq = 0;
        ip;
        debugflag = false;
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
        reset_stm = 9;
        stm_ready = 10;
        hold_input = 11;


        stop_sound = 0;
        start_tone = 1;
        start_sweep = 2;
        start_noise = 3;
        start_ripple = 4;
        start_wav = 5;
        start_ba = 6;
    end

    methods
        function this = rz6_unircx_tasklist.m(rz6_module)
            this.module = rz6_module;
            this.tag_name = 'wherethedatagoes';
        end
        
        function tq_i32=get(this)
            tq_i32 = int32(this.tq(1:this.nq,:));     
        end
        
 
        function add_task(this, delaytime, task, varargin)
           desc = this.parse_command(delaytime, task, varargin{:});
           this.append_task(desc);
        end

        function debug(this, val)
           if (val)
              this.debugflag = true;
           else
              this.debugflag = false;
           end
        end
    end


    methods (Access=protected)
        function append_task(this, desc)
           this.nq = this.nq+1;
           this.tq(this.nq,:) = struct2array(desc);
        end

        function desc = parse_command(this, varargin)
            funcName='rz6_unircx_tasklist.m/add_task';
            validCommands = { 'WaitForTrigger', 'SoundA','SoundB','Mux','Signaling',...
               'SoundMov','Daq','SetDIO','TrigOut','Reset','Ready','HoldInp'};

            validCommand = @(x) any(validatestring(x,validCommands));
            validDelayTime = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});

            p = inputParser;
            p.FunctionName = funcName;
            p.addRequired('DelayTime', validDelayTime);
            p.addRequired('Command', validCommand);
            p.parse(varargin{1:2});

            switch lower(p.Results.Command)
                case 'waitfortrigger'
                    desc = this.parse_waitfortrigger(p,varargin{:});

                case 'sounda'
                    desc = this.parse_sounda(p,varargin{:});

                case 'soundb'
                    desc = this.parse_soundb(p,varargin{:});

                case 'mux'
                    desc = this.parse_mux(p,varargin{:});

                case 'signaling'
                    desc = this.parse_signaling_byte(p,varargin{:});

                case 'setsignalingbyte'
                    desc = this.parse_setsignalingbyte(p,varargin{:});

                case 'soundmov'
                    desc = this.parse_soundmov(p,varargin{:});

                case 'daq'
                    desc = this.parse_daq(p,varargin{:});

                case 'setdio'
                    desc = this.parse_setdio(p,varargin{:});

                case 'trigout'
                    desc = this.parse_trigout(p,varargin{:});

                case 'reset'
                    desc = this.parse_reset(p,varargin{:});

                case 'ready'
                    desc = this.parse_ready(p,varargin{:});

                case 'holdinp'
                    desc = this.parse_holdinp(p,varargin{:});

            otherwise
                error('unknown task, this is a bug');
            end

            desc.DelayTime = p.Results.DelayTime;

            if this.debugflag
               disp(['---',mfilename,' debug output---']);
               disp('INPUT:');
               disp(p.Results);
               disp('');
               disp('RESULT:');
               disp(desc);
               disp('---end---');
            end

            if any(isnan(struct2array(desc)))
               error('missing arguments');
            end
        end

        function desc = newdesc(~)
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

        function desc = parse_waitfortrigger(this,p,varargin)
            desc = this.newdesc();
            desc.TaskType = this.wait_for_trigger;

            validExtTrig = @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative','<',8});
            
            p.addRequired('Input', @(x) validInput(x));
            p.addOptional('ExternalTrigger', 0, @(x) validExtTrig(x));
            p.parse(varargin{:});
            
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

        function desc = parse_sounda(this,p,varargin)
            desc = this.parse_sound(1,p,varargin{:});
        end

        function desc = parse_soundb(this,p,varargin)
            desc = this.parse_sound(2,p,varargin{:});
        end

        function desc = parse_sound(this,channel,p,varargin)
            desc = this.newdesc();
            if channel==1
               desc.TaskType = this.start_stop_sound_a;
            elseif channel==2
               desc.TaskType = this.start_stop_sound_b;
            else
               error('expected channel to be 1 or 2, this is a bug');
            end

            validFreq = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
            validAttenuation = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
            validPosNum = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
            validStartPhase = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-180,'<=',180});
            expectedSounds = { 'None','Stop','Tone','Sweep','Noise','Ripple','Wav','B=A'};

            p.addRequired('SoundType', @(x) any(validatestring(x, expectedSounds)));
            p.parse(varargin{1:3});

            % Expand partially matched SoundType strings, and convert to lowercase
            SoundType=lower(validatestring(p.Results.SoundType, expectedSounds));

            switch SoundType
            case 'stop'
               desc.SoundType = 0;
               p.parse(varargin{:});

            case 'tone'
               desc.SoundType = 1;
               p.addRequired('fCenter', validFreq);
               p.addRequired('fModulation', validFreq);
               p.addRequired('bwModulation', validFreq);
               p.addRequired('Attenuation', validAttenuation);
               p.parse(varargin{:});
               desc.Par1 = p.Results.fCenter;
               desc.Par2 = p.Results.fModulation; % f_mod
               desc.Par3 = p.Results.bwModulation; % mod_bw
               desc.Par4 = p.Results.Attenuation; % att

            case 'sweep'
               desc.SoundType = 2;
               p.addRequired('fStart', validFreq);
               p.addRequired('nOctaves', validPosNum);
               p.addRequired('nSweeps', validPosNum);
               p.addRequired('Attenuation', validAttenuation);
               p.parse(varargin{:});
               desc.Par1 = p.Results.fStart;
               desc.Par2 = p.Results.nOctaves;
               desc.Par3 = p.Results.nSweeps; % sweeps/1000sec
               desc.Par4 = p.Results.Attenuation;

            case 'noise'
               desc.SoundType = 3;
               p.addRequired('fLowPass', validFreq);
               p.addRequired('fHighPass', validFreq);
               p.addRequired('Attenuation', validAttenuation);
               p.parse(varargin{:});
               desc.Par1 = p.Results.fLowPass;
               desc.Par2 = p.Results.fHighPass;
               desc.Par4 = p.Results.Attenuation;

            case 'ripple'
               desc.SoundType = 4;
               p.addRequired('fStart', validFreq);
               p.addRequired('modInTime', validPosNum);
               p.addRequired('modInFreq', validPosNum);
               p.addRequired('Attenuation', validAttenuation);
               p.parse(varargin{:});
               desc.Par1 = p.Results.fStart;
               desc.Par2 = p.Results.modInTime; % Hz
               desc.Par3 = p.Results.modInFreq; % Phase/Octave
               desc.Par4 = p.Results.Attenuation;

            case 'wav'
               desc.SoundType = 5;
               p.addRequired('wavIndex', validPosNum);
               p.addRequired('Attenuation', validAttenuation);
               p.parse(varargin{:});
               desc.Par1 = p.Results.wavIndex;
               desc.Par4 = p.Results.Attenuation;

            case 'b=a'
               desc.SoundType = 6;
               expectedMovTypes = { 'None','Linear','Sine' };
               p.addRequired('movType', @(x) any(validatestring(x, expectedMovTypes)));
               p.parse(varargin{1:4});
               switch lower(p.Results.movType)
               case 'none'
                  p.addRequired('Attenuation', validAttenuation);
                  p.parse(varargin{:});
                  desc.Par1 = 0;
                  desc.Par4 = p.Results.Attenuation;
               case 'sine'
                  p.addRequired('tPeriod', validPosNum);
                  p.addRequired('startPhase', validStartPhase);
                  p.addRequired('Attenuation', validAttenuation);
                  p.parse(varargin{:});
                  desc.Par1 = 1;
                  desc.Par2 = p.Results.tPeriod; % msec
                  desc.Par3 = p.Results.startPhase; 
                  desc.Par4 = p.Results.Attenuation;
               case 'linear'
                  p.addRequired('tPeriod', validPosNum);
                  p.addRequired('startPhase', validStartPhase);
                  p.addRequired('Attenuation', validAttenuation);
                  p.parse(varargin{:});
                  desc.Par1 = 2;
                  desc.Par2 = p.Results.tPeriod; % msec
                  desc.Par3 = p.Results.startPhase; 
                  desc.Par4 = p.Results.Attenuation;
               otherwise
                  error('invalid mov type, this is a bug');
               end

            otherwise
               error('unknown sound type, this is a bug');
            end

        end

        function desc = parse_mux(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.set_mux;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           p.addRequired('MuxByte', @(x) validByte(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.MuxByte;
        end

        function desc = parse_signaling_byte(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.set_signaling_byte;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           p.addRequired('SignalingByte', @(x) validByte(x));
           p.parse(varargin{:});
           desc.Par1 = p.Results.SignalingByte;
        end

        function desc = parse_soundmov(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.start_stop_moving_sounds;

           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));
           validNumSpeakers = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','>',1});
           validPeriod = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           validStartPhase = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>=',-180, '<=',180});

           p.addRequired('StartStop', validStartStop);
           p.parse(varargin{1:3});
           switch lower(p.Results.StartStop)
           case 'stop'
               p.parse(varargin{:});
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
               p.addRequired('NumSpeakers', validNumSpeakers);
               p.addRequired('Period', validPeriod);
               p.addRequired('StartPhase', validStartPhase);
               p.parse(varargin{:});
               desc.Par2 = p.Results.NumSpeakers;
               desc.Par3 = p.Results.Period;
               desc.Par4 = p.Results.StartPhase;
           otherwise
               error('unexcpected error, this is a bug');
           end
        end

        function desc = parse_daq(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.start_stop_daq;

           validInputSelByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           validDivisor = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','positive'});
           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));

           p.addRequired('StartStop', validStartStop);
           p.parse(varargin{1:3});
           switch lower(p.Results.StartStop)
           case 'stop'
               p.parse(varargin{:});
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
               p.addRequired('InputSelByte', validInputSelByte);
               p.addOptional('Divisor', 1, validDivisor);
               p.parse(varargin{:});
               desc.Par2 = p.Results.InputSelByte;
               desc.Par3 = p.Results.Divisor;
           otherwise
               error('unexcpected error, this is a bug');
           end
        end

        function desc = parse_setdio(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.set_digital_out;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});

           p.addRequired('OutputByte', validByte);
           p.parse(varargin{:});
           desc.Par1 = p.Results.OutputByte;
        end

        function desc = parse_trigout(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.output_trigger;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           validDTInterval = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>',100});


           p.addRequired('OutputByte', validByte);
           p.addOptional('DoubleTriggerInterval',0, validDTInterval);
           p.parse(varargin{:});
           desc.Par1 = p.Results.OutputByte;
           desc.Par2 = p.Results.DoubleTriggerInterval;
        end

        function desc = parse_reset(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.reset_stm;
           p.parse(varargin{:});
        end

        function desc = parse_ready(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.stm_ready;
           p.parse(varargin{:});
        end

        function desc = parse_holdinp(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.hold_input;

           validInputMask = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative', '<',256});
           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));

           p.addRequired('StartStop', validStartStop);
           p.parse(varargin{1:3});
           switch lower(p.Results.StartStop)
           case 'stop'
               p.parse(varargin{:});
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
               p.addRequired('InputMask', validInputMask);
               p.parse(varargin{:});
               desc.Par2 = p.Results.InputMask;
           otherwise
               error('unexcpected error, this is a bug');
           end

        end


    end

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
