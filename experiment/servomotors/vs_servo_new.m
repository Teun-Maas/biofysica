classdef vs_servo_new < panaservo

% 2018-01-02 GW: fixed maxint problem in convert_profile() function.
%            Now converting position values to int32 instead of int16.
% 


methods (Access=protected)
    
    function set_enable(this, value)
        this.plc.IEC_write(this.varmap.Mch_Enable_Servo_1, value==1);
        this.plc.IEC_write(this.varmap.Mch_Enable_Servo_2, value==1);
    end
    
end

methods
    function this = vs_servo_new
        this@panaservo('192.168.1.10', 'vestibular_PLC_global_variables.xlsx');
    end

    function r = is_enabled(this)
        r = ...
           this.plc.IEC_read(this.varmap.Mch_Enable_Servo_1) & ...
           this.plc.IEC_read(this.varmap.Mch_Enable_Servo_2);
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
            'X_Emergency'
            'X_Key_Switch'
            'X_SW4_S1_Home'
            'X_SW3_S2_Home'
            'X_SW1_S2_Limit_CW'
            'X_SW2_S2_Limit_CCW'
            'X_SW5_Port_Closed_Left'
            'X_SW7_Port_Closed_Right'
            'X_SD3A1_No_Alarm_1'
            'X_SD3A1_No_Alarm_1_Rear'
            'Mch_Enable_Servo_1'
            'Mch_Enable_Servo_2'
            'Servo_1_Enabled'
            'Servo_2_Enabled'
            'Run_Axis1_Prg_A'
            'Run_Axis1_Prg_B'
            'Run_Axis1_Prg_C'
            'Run_Axis2_Prg_A'
            'Run_Axis2_Prg_B'
            'Run_Axis2_Prg_C'
            'ContinueWith_T1A'
            'ContinueWith_T1B'
            'ContinueWith_T1C'
            'ContinueWith_T2A'
            'ContinueWith_T2B'
            'ContinueWith_T2C'
            'cmd_Run_Program_A'
            'cmd_Run_Program_B' 
            'cmd_Run_Program_C'
            'PV_Table_1A_Index'
            'PV_Table_1B_Index'
            'PV_Table_1C_Index'
            'PV_Table_2A_Index'
            'PV_Table_2B_Index'
            'PV_Table_2C_Index'
            %'PV_Position_1_Pls'
            %'PV_Position_2_Pls'
        };
        for i = 1:length(status_vars)
            v=status_vars{i};
            vd=this.varmap.lookup(v);
            val=this.plc.IEC_read(vd);
            %fprintf('%c[K%s:\t%d\n', 27, v, val);
            fprintf('%s:\t%d\n', v, val);

        end
    end

    function result = write_profile(this, axis1, axis2)            
       p1 = this.convert_profile(axis1);
       p2 = this.convert_profile(axis2);
       this.plc.IEC_write(this.varmap.Table_1A, p1);
       this.plc.IEC_write(this.varmap.Table_2A, p2);
       result = 0;
    end

    function [ax1, ax2] = read_profile_sv(this)
       ax1=this.plc.IEC_read(this.varmap.Table_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.Table_2A, 0, 2000);
       ax1=cumsum(double(ax1))/10.0;
       ax2=cumsum(double(ax2))/10.0;
    end

    function [ax1, ax2] = read_profile_pv(this)
       ax1=this.plc.IEC_read(this.varmap.PV_Position_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.PV_Position_2A, 0, 2000);
       ax1=double(ax1)/10.0;
       ax2=double(ax2)/10.0;
    end

end  % methods

end  % classdef
