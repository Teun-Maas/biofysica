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
        
        function append_task(this, desc)
           this.nq = this.nq+1;
           this.tq(this.nq,:) = desc;
        end

    end


    methods (Access=protected)
        function desc = parse_command(this,command,varargin)
            switch lower(command)
                case 'waitfortrigger'
                    desc = this.parse_waitfortrigger(varargin);

                case 'startsound'
                    desc = this.parse_startsound(varargin);

                case 'stopsound'
                    desc = this.parse_stopsound(varargin);

                case 'setmux'
                    desc = this.parse_setmux(varargin);

                case 'resetmux'
                    desc = this.parse_resetmux(varargin);

                case 'setsignalingbyte'
                    desc = this.parse_setsignalingbyte(varargin);

                case 'startmovingsound'
                    desc = this.parse_startmovingsound(varargin);

                case 'stopmovingsound'
                    desc = this.parse_stopmovingsound(varargin);

                case 'startdaq'
                    desc = this.parse_startdaq(varargin);

                case 'stopdaq'
                    desc = this.parse_stopdaq(varargin);

                case 'setdigitalout'
                    desc = this.parse_setdigitalout(varargin);

                case 'outputtrigger'
                    desc = this.parse_outputtrigger(varargin);

            otherwise
                error('unknown task');
            end
        end

        function desc = newdesc(this)
            desc = struct(
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
            validDelayTime = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x < 2^31);
            validExternalTrigger = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x < 8);
            p.addRequired('DelayTime', validDelayTime);
            p.addRequired('Input', @(x) any(validatestring(x,expectedInputs)));
            p.addOptional('ExternalTrigger', validByte);
            desc(4) = p.Results.input;

            function validInput(x)
               stringInputs = { 'ZBusB', 'External', 'Soft1' };
               if isstring(x)
            end
        end

        function desc = parse_startsound(this,varargin)
        end

        function desc = parse_stopsound(this,varargin)
        end

        function desc = parse_setmux(this,varargin)
        end

        function desc = parse_resetmux(this,varargin)
        end

        function desc = parse_setsignalingbyte(this,varargin)
        end

        function desc = parse_startmovingsound(this,varargin)
        end

        function desc = parse_stopmovingsound(this,varargin)
        end

        function desc = parse_startdaq(this,varargin)
        end

        function desc = parse_stopdaq(this,varargin)
        end

        function desc = parse_setdigitalout(this,varargin)
        end

        function desc = parse_outputtrigger(this,varargin)
        end

        %function desc = parse_(this,varargin)
        %end

        %function desc = parse_(this,varargin)
        %end


    end

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
