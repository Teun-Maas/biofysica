classdef vs_servo < handle

% 2018-01-02 GW: fixed maxint problem in convert_profile() function.
%            Now converting position values to int32 instead of int16.
% 

properties
    varmap;
    plc;

end

methods
    function this = vs_servo
        this.plc = m2c_plc('192.168.1.10');
        this.varmap = m2c_memory_map('vestibular_PLC_global_variables.xlsx');
    end

    function delete(this)
        delete(this.varmap);
        delete(this.plc);
    end

    function set_enable(this, value)
        this.plc.IEC_write(this.varmap.Mch_Enable_Servo_1, value==1);
        this.plc.IEC_write(this.varmap.Mch_Enable_Servo_2, value==1);
    end
    
    function enable(this)
        if ~this.is_enabled()
            warndlg(sprintf('Enable servo drives, then\npress OK to continue'),...
                'Warning','modal');
        end
    end
    
    function disable(this)
        this.set_enable(0);
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

    function status = print_status(this)
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

    function idata = convert_profile(this, data, limit_r_vel, limit_r_pos)
        maxlen = 2000;

        data = reshape(data, [], 1); % make it a column vector
        len = length(data);
        if len > maxlen
           warning('profile too long, truncated to %d samples', maxlen);
           data = data(1:maxlen);
        end
        npad = maxlen-len;
        % convert to int32 data in 1/10 degrees, padded with zeros
        data = int32(round(10*data));
        r_pos = 10*range(data);  % r_pos in degrees
        if  nargin >=4
           if (r_pos(1) < limit_r_pos(1)) | (r_pos(2) > limit_r_pos(2))
              warning('profile position out of range');
           end
        end
% FIXME?
% ??? pad with one leading zero before d/dt? to keep buffer length == 2000
% padding is fine when running profile from only one buffer.

        idata = diff(data);
        idata = [0; idata];

        if nargin >= 3
            r_vel = max(abs(range(idata)));  % r_vel in degrees/second
            if (r_vel > limit_r_vel)
               warning('profile speed out of range');
            end
        end
        idata = int16([idata; zeros(npad, 1)]);
    end
    
    function result = clear_profile(this);
       result = this.write_profile([],[]);
    end

    function result = write_profile(this, axis1, axis2);     
        
       p1 = this.convert_profile(axis1);
       p2 = this.convert_profile(axis2);
       
       this.plc.IEC_write(this.varmap.Table_1A, p1);
       this.plc.IEC_write(this.varmap.Table_2A, p2);
       result = 0;
    end

    function [ax1, ax2] = read_profile_sv(this);
       ax1=this.plc.IEC_read(this.varmap.Table_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.Table_2A, 0, 2000);
       ax1=cumsum(double(ax1))/10.0;
       ax2=cumsum(double(ax2))/10.0;
    end

    function [ax1, ax2] = read_profile_pv(this);
       ax1=this.plc.IEC_read(this.varmap.PV_Position_1A, 0, 2000);
       ax2=this.plc.IEC_read(this.varmap.PV_Position_2A, 0, 2000);
       ax1=double(ax1)/10.0;
       ax2=double(ax2)/10.0;
    end

end  % methods

end  % classdef
