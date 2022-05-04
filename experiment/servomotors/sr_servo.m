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
        
        function idata = convert_profile(~, data, limit_r_vel, limit_r_pos)
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
                if (r_pos(1) < limit_r_pos(1)) || (r_pos(2) > limit_r_pos(2))
                    warning('profile position out of range');
                end
            end
            % FIXME?
            % ??? pad with one leading zero before d/dt? to keep buffer length == 2000
            % padding is fine when running profile from only one buffer.
            
            % GW: Workaround 2021-11-01 / FIXED 2022-05-02
                         idata = diff(data);
                         idata = [0; idata];
            %idata=data;
            %
            if nargin >= 3
                r_vel = max(abs(range(idata)));  % r_vel in degrees/second
                if (r_vel > limit_r_vel)
                    warning('profile speed out of range');
                end
            end
            idata = int16([idata; zeros(npad, 1)]);
        end
        
        function result = write_profile(this, main, speaker, chair)
            p1 = this.convert_profile(main);
            p2 = this.convert_profile(chair);
            p3 = this.convert_profile(speaker);
          
            this.plc.IEC_write(this.varmap.Profile_1A, p1);
            this.plc.IEC_write(this.varmap.Profile_2A, p2);
            this.plc.IEC_write(this.varmap.Profile_3A, p3);
            result = 0;
        end
        
        function [main, speaker, chair] = read_profile_sv(this)
            main=this.plc.IEC_read(this.varmap.Profile_1A, 0, 2000);
            chair=this.plc.IEC_read(this.varmap.Profile_2A, 0, 2000);
            speaker=this.plc.IEC_read(this.varmap.Profile_3A, 0, 2000);
            % GW: workaround 2021-11-01 / FIXED 2022-05-02
            main=cumsum(double(main))/10.0;
            chair=cumsum(double(chair))/10.0;
            speaker=cumsum(double(speaker))/10.0;
        end
        
        function [main, speaker, chair] = read_profile_pv(this)
            main=this.plc.IEC_read(this.varmap.PV_Position_1A, 0, 2000);
            chair=this.plc.IEC_read(this.varmap.PV_Position_2A, 0, 2000);
            speaker=this.plc.IEC_read(this.varmap.PV_Position_3A, 0, 2000);
            % GW FIXME
%             main=double(main)/115;
%             chair=double(chair)/115;
%             speaker=double(speaker)/10;
            main=double(main)/10;
            chair=double(chair)/10;
            speaker=double(speaker)/10;
        end
        
    end  % methods
    
end  % classdef
