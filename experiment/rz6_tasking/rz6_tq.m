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
        set_reset_mux = 3;
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
        
 
        function add_task(this, task, varargin)
           desc = this.parse_command(task, varargin{:});
           this.append_task(desc);
        end
    end


    methods (Access=protected)
        function append_task(this, desc)
           this.nq = this.nq+1;
           this.tq(this.nq,:) = struct2array(desc);
        end

        function desc = parse_command(this,command, delaytime, varargin)
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
        end

        function desc = parse_soundb(this,varargin)
        end

        function desc = parse_mux(this,varargin)
        end

        function desc = parse_signaling(this,varargin)
        end

        function desc = parse_setsignalingbyte(this,varargin)
        end

        function desc = parse_soundmov(this,varargin)
        end

        function desc = parse_daq(this,varargin)
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
