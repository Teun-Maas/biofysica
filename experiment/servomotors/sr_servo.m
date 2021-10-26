classdef sr_servo < panaservo


methods (Access=protected)
    
    function set_enable(this, value)
        this.plc.IEC_write(this.varmap.mch_Enable_Servo_1, value==1);
        this.plc.IEC_write(this.varmap.mch_Enable_Servo_2, value==1);
        this.plc.IEC_write(this.varmap.mch_Enable_Servo_3, value==1);
    end
    
end

methods
    function this = sr_servo
        this@panaservo('192.168.1.5', 'speakerrobot_PLC_global_variables.xlsx');
    end

    function r = is_enabled(this)
        r = ...
           this.plc.IEC_read(this.varmap.mch_Enable_Servo_1) & ...
           this.plc.IEC_read(this.varmap.mch_Enable_Servo_2) & ...
           this.plc.IEC_read(this.varmap.mch_Enable_Servo_3);
    end
    
    function start(this)
        this.plc.IEC_write(this.varmap.cmd_Run_Program_A, 1);
    end

    function stop(this)
        this.plc.IEC_write(this.varmap.cmd_Run_Program_A, 0);
    end

    function print_status(this)
        status_vars = {
            'X_S1_Ready'
            'X_S1_Alarm'
            'X_S2_Ready'
            'X_S2_Alarm'
            'X_S3_Ready'
            'X_S3_Alarm'
            'X_No_Emergency'
            'X_Key_Unlock'
            'mch_Enable_Servo_1'
            'mch_Enable_Servo_2'
            'mch_Enable_Servo_3'
            'Y_S1_Servo_ON'
            'Y_S2_Servo_ON'
            'Y_S3_Servo_ON'
            'Run_Axis1_Prg_A'
            'Run_Axis1_Prg_B'
            'Run_Axis1_Prg_C'
            'Run_Axis2_Prg_A'
            'Run_Axis2_Prg_B'
            'Run_Axis2_Prg_C'
            'Run_Axis3_Prg_A'
            'Run_Axis3_Prg_B'
            'Run_Axis3_Prg_C'
            'cmd_Run_Program_A'
            'cmd_Run_Program_B' 
            'cmd_Run_Program_C'
        };
        for i = 1:length(status_vars)
            v=status_vars{i};
            vd=this.varmap.lookup(v);
            val=this.plc.IEC_read(vd);
            %fprintf('%c[K%s:\t%d\n', 27, v, val);
            fprintf('%s:\t%d\n', v, val);

        end
    end

    function result = write_profile(this, axis1, axis2, axis3)            
       p1 = this.convert_profile(axis1);
       p2 = this.convert_profile(axis2);
       p3 = this.convert_profile(axis3);
       this.plc.IEC_write(this.varmap.Profile_1A, p1);
       this.plc.IEC_write(this.varmap.Profile_2A, p2);
       this.plc.IEC_write(this.varmap.Profile_3A, p3);
       result = 0;
    end

    function [ax1, ax2, ax3] = read_profile_sv(this)
       ax1=this.plc.IEC_read(this.varmap.Profile_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.Profile_2A, 0, 2000);
       ax3=this.plc.IEC_read(this.varmap.Profile_3A, 0, 2000);
       ax1=cumsum(double(ax1))/10.0;
       ax2=cumsum(double(ax2))/10.0;
       ax3=cumsum(double(ax3))/10.0;
    end

    function [ax1, ax2, ax3] = read_profile_pv(this)
       ax1=this.plc.IEC_read(this.varmap.PV_Position_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.PV_Position_2A, 0, 2000);
       ax3=this.plc.IEC_read(this.varmap.PV_Position_3A, 0, 2000);
       ax1=double(ax1)/10.0;
       ax2=double(ax2)/10.0;
       ax3=double(ax3)/10.0;
    end

end  % methods

end  % classdef
