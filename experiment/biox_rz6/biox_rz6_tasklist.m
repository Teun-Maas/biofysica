classdef biox_rz6_tasklist < handle
    properties (Access=protected)
        tq = [];
        nq = 0;
        ip;
        debugflag = false;
    end

    properties (Constant)
        desclen = 7;

        task_waitfortrigger = 0;
        task_sound_a = 1;
        task_sound_b = 2;
        task_mux = 3;
        task_holdinput = 4; 
        task_soundmov = 5;
        task_daq = 6;
        task_dioout = 7;
        task_trgout = 8;
        task_reset = 9;
        tsk_ready = 10;
        task_multiconfiga = 11;
        task_multiconfigb = 12;
        task_atta = 13;
        task_attb = 14;

        soundtype_stop = 0;
        soundtype_tone = 1;
        soundtype_sweep = 2;
        soundtype_noise = 3;
        soundtype_ripple = 4;
        soundtype_wav = 5;
        soundtype_multitone = 6;
        soundtype_ba = 7;
    end

    methods
        function this = biox_rz6_tasklist 
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
            funcName='biox_rz6_tasklist.m/add_task';
            validCommands = { 'WaitForTrigger', 'SoundA','SoundB','Mux','HoldInput',...
               'SoundMov','Daq','SetDIO','TrigOut','Reset','Ready',...
               'MultiConfigA','MultiConfigB','AttA','AttB'};

            validCommand = @(x) any(validatestring(x,validCommands));
            validDelayTime = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});

            p = biox_inputParser;
            p.FunctionName = funcName;
            p.addRequired('DelayTime', validDelayTime);
            p.addRequired('Command', validCommand);
            p.parse(varargin{1:2});

            % Expand partially matched commands 
            command = lower(validatestring(p.Results.Command,validCommands));
            switch command
                case 'waitfortrigger'
                    desc = this.parse_waitfortrigger(p,varargin{:});

                case 'sounda'
                    desc = this.parse_sounda(p,varargin{:});

                case 'soundb'
                    desc = this.parse_soundb(p,varargin{:});

                case 'mux'
                    desc = this.parse_mux(p,varargin{:});

                case 'holdinput'
                    desc = this.parse_holdinput(p,varargin{:});

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

                case 'task_multiconfiga'
                    desc = this.parse_multiconfiga(p,varargin{:});

                case 'task_multiconfigb'
                    desc = this.parse_multiconfigb(p,varargin{:});

                case 'atta'
                    desc = this.parse_atta(p,varargin{:});

                case 'attb'
                    desc = this.parse_attb(p,varargin{:});

            otherwise
                error('unknown task: %s', p.Results.Command);
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
            desc.TaskType = this.task_waitfortrigger;

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
               desc.TaskType = this.task_sound_a;
            elseif channel==2
               desc.TaskType = this.task_sound_b;
            else
               error('expected channel to be 1 or 2, this is a bug');
            end

            validFreq = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
            validPosNum = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});
            valid01 = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1});
            validStartPhase = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-180,'<=',180});
            expectedSounds = { 'None','Stop','Tone','Sweep','Noise','Ripple','Wav','B=A'};

            p.addRequired('SoundType', @(x) any(validatestring(x, expectedSounds)));
            p.parse(varargin{1:3});

            % Expand partially matched SoundType strings, and convert to lowercase
            SoundType=lower(validatestring(p.Results.SoundType, expectedSounds));

            switch SoundType
            case 'stop'
               desc.SoundType = this.soundtype_stop;
               p.parse(varargin{:});

            case 'tone'
               desc.SoundType = this.soundtype_tone;
               p.addRequired('ToneFreq', validFreq);
               p.addRequired('ModFreq', validFreq);
               p.addRequired('ModBW', validFreq);
               p.parse(varargin{:});
               desc.Par1 = p.Results.ToneFreq;
               desc.Par2 = p.Results.ModFreq;
               desc.Par3 = p.Results.ModBW;

            case 'sweep'
               desc.SoundType = this.soundtype_sweep;
               p.addRequired('StartFreq', validFreq);
%TODO change in manual, NrOctaves vs. no_octaves
               p.addRequired('NrOctaves', validPosNum);
               p.addRequired('Period', validPosNum);
               p.parse(varargin{:});
               desc.Par1 = p.Results.StartFreq;
               desc.Par2 = p.Results.NrOctaves;
               desc.Par3 = p.Results.Period; % Period in msec

            case 'noise'
               desc.SoundType = this.soundtype_noise;
               p.addRequired('HpFreq', validFreq);
               p.addRequired('LpFreq', validFreq);
               p.parse(varargin{:});
               desc.Par1 = p.Results.HpFreq;
               desc.Par2 = p.Results.LpFreq;

            case 'ripple'
               desc.SoundType = this.soundtype_ripple;
               p.addRequired('StartFreq', validFreq);
               p.addRequired('ModInTime', validPosNum);
               p.addRequired('ModInFreq', validPosNum);
               p.parse(varargin{:});
               desc.Par1 = p.Results.StartFreq;
               desc.Par2 = p.Results.ModInTime; % Hz
               desc.Par3 = p.Results.ModInFreq; % Phase/Octave

            case 'wav'
               desc.SoundType = this.soundtype_wav;
               p.addRequired('Reset', valid01);
               p.parse(varargin{:});
               desc.Par1 = p.Results.Reset;

            case 'multitone'
               desc.SoundType = this.soundtype.multitone;
               p.parse(varargin{:});

            case 'b=a'
               assert(desc.TaskType == this.task_sound_b,...
                  'SoundType ''b=a'' is only valid for taskType ''SoundB''');
               desc.SoundType = this.soundtype_ba;
               expectedMovTypes = { 'None','Linear','Sine' };
               p.addRequired('movType', @(x) any(validatestring(x, expectedMovTypes)));
               p.parse(varargin{1:4});
               movType=lower(validatestring(p.Results.movType,expectedMovTypes));
               switch lower(movType)
               case 'none'
                  p.parse(varargin{:});
                  desc.Par1 = 0;
               case 'sine'
                  p.addRequired('Period', validPosNum);
                  p.addRequired('Phase', validStartPhase);
                  p.parse(varargin{:});
                  desc.Par1 = 1;
                  desc.Par2 = p.Results.Period; % msec
                  desc.Par3 = p.Results.Phase; 
               case 'linear'
                  p.addRequired('Period', validPosNum);
                  p.addRequired('Phase', validStartPhase);
                  p.parse(varargin{:});
                  desc.Par1 = 2;
                  desc.Par2 = p.Results.Period; % msec
                  desc.Par3 = p.Results.Phase; 
               otherwise
                  error('invalid mov type, this is a bug');
               end

            otherwise
               error('unknown sound type, this is a bug');
            end

        end

        function desc = parse_mux(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_mux;

           validChan = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',15});
           validDevice = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',3});
           % validSet = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',1});
           % validReset = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',1});
           expectedActions = { 'Set','Reset' };
           p.addRequired('Device', @(x) validDevice(x));
           p.addRequired('Action', @(x) any(validatestring(x, expectedActions)));
           p.addOptional('Channel', [], @(x) validChan(x));

           % p.addRequired('Set', @(x) validSet(x));
           % p.addRequired('Reset', @(x) validReset(x));
           p.parse(varargin{:});

           if strcmpi(p.Results.Action,'Set')
               if isempty(p.Results.Channel)
                   error('''Set'' requires a the value of ''Channel'' to be specified.');
               else
                   Channel = p.Results.Channel;
               end
               ASet=1; AReset=0;
           elseif strcmpi(p.Results.Action,'Reset')
               Channel = 0;
               ASet=0; AReset=1;
           else
               error('Unexpected or missing ''Action'', this is a bug');
           end   
           desc.Par1 = Channel + 16*p.Results.Device + 64*ASet + 128*AReset;
        end

        function desc = parse_holdinput(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_holdinput;

           validInputMask = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative', '<',256});

           p.addRequired('InputMask', validInputMask);
           p.parse(varargin{:});
           desc.Par1 = p.Results.InputMask;
        end

        function desc = parse_soundmov(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_soundmov;

           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));
           validNumSpeakers = @(x) validateattributes(x, ...
              {'numeric'}, {'odd','scalar','>=',3,'<=',21});
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
           desc.TaskType = this.task_daq;

           validChannelSelection = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           validDivisor = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','positive'});
           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));

           p.addRequired('ChannelSelection', validChannelSelection);
           p.addRequired('Divisor', validDivisor);
           p.addRequired('StartStop', validStartStop);

% % TODO new version, always have three parameters
           p.parse(varargin{:});
           switch lower(p.Results.StartStop)
           case 'stop'
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
           otherwise
               error('unexcpected error, this is a bug');
           end
           desc.Par2 = p.Results.ChannelSelection;
           desc.Par3 = p.Results.Divisor;

% % TODO This is the previous version. Delete if above is correct           
% %            p.addRequired('StartStop', validStartStop);
% %            p.parse(varargin{1:3});
% %            switch lower(p.Results.StartStop)
% %            case 'stop'
% %                p.parse(varargin{:});
% %                desc.Par1 = 0;
% %            case 'start'
% %                desc.Par1 = 1;
% %                p.addRequired('ChannelSelection', validChannelSelection);
% %                p.addOptional('Divisor', 1, validDivisor);
% %                p.parse(varargin{:});
% %                desc.Par2 = p.Results.ChannelSelection;
% %                desc.Par3 = p.Results.Divisor;
% %            otherwise
% %                error('unexcpected error, this is a bug');
% %            end
        end

        function desc = parse_setdio(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_dioout;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});

           p.addRequired('OutputByte', validByte);
           p.parse(varargin{:});
           desc.Par1 = p.Results.OutputByte;
        end

        function desc = parse_trigout(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_trgout;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           validDTInterval = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>',100});


           p.addRequired('OutputByte', validByte);
%TODO change in manual: DoubleDelay vs DelayTime
           p.addOptional('DoubleDelay',0, validDTInterval);
           p.parse(varargin{:});
           desc.Par1 = p.Results.OutputByte;
           desc.Par2 = p.Results.DoubleDelay;
        end

        function desc = parse_reset(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_reset;
           p.parse(varargin{:});
        end

        function desc = parse_ready(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.tsk_ready;
           p.parse(varargin{:});
        end

        function desc = parse_multiconfig(this,tasktype,p,varargin)
           validIndex = @(x) validateattributes(x,{'numeric'},{'scalar','positive','<=',4});
           validFrequency = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
           validPhase = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-180,'<=',180});
%TODO check manual if 0 <= amplitude <=1 is correct
           validAmplitude = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1});

           desc = this.newdesc();
           desc.TaskType = tasktype;
           p.addRequired('Index', validIndex);
           p.addRequired('Frequency', validFrequency);
           p.addRequired('Phase', validPhase);
           p.addRequired('Amplitude', validAmplitude);
           p.parse(varargin{:});
           desc.Par1 = p.Results.Index;
           desc.Par2 = p.Results.Frequency;
           desc.Par3 = p.Results.Phase;
           desc.Par4 = p.Results.Amplitude;
        end

        function desc = parse_multiconfiga(this,p,varargin)
            desc = this.parse_multiconfig(this.task_multiconfiga, p, varargin{:});
        end

        function desc = parse_multiconfig_b(this,p,varargin)
            desc = this.parse_multiconfig(this.task_multiconfigb, p, varargin{:});
        end

        function desc = parse_att(this,tasktype,p,varargin)
           validAttenuation = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',80});
           desc = this.newdesc();
           desc.TaskType = tasktype;
           p.addRequired('Attenuation', validAttenuation);
           p.parse(varargin{:});
           desc.Par1 = p.Results.Attenuation;
        end

        function desc = parse_atta(this,p,varargin)
            desc = this.parse_att(this.task_atta,p,varargin{:});
        end

        function desc = parse_attb(this,p,varargin)
            desc = this.parse_att(this.task_attb,p,varargin{:});
        end

    end

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
