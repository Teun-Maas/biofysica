classdef rz6_tq < handle
    properties (Access=protected)
        module;
        tag_name;
        tq = [];
        nq = 0;
    end

    properties (Constant)
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

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
