classdef rz6_task_queue < handle
    properties (Access=protected)
        module;
        tag_name;
        tq = [];
        nq = 0;
    end

    methods
        function this = rz6_task_queue(rz6_module)
            this.module = rz6_module;
            this.tag_name = 'wherethedatagoes';
        end
        
        function upload(this)
            tq_i32 = int32(this.tq(1:this.nq,:))';     
            nOS = 0;
            %TODO
            this.module.WriteTagVEX(this.tag_name, nOS, 'I32', tq_i32);
        end
        
        function append_task(this, task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
           this.nq = this.nq+1;
           this.tq(this.nq,:) = desc;
        end
        
    end
end
