classdef m2c_plc < handle

properties (Access=protected)

   handle;
   wd_timestamp;
   wd_timer;
   wd_keepalive_interval;
end

methods

    function this = m2c_plc(ipAddress)
        this.handle = m2c_new_tcp(ipAddress, 9094);
        this.wd_timestamp = tic;
        this.wd_keepalive_interval = 30;
        this.wd_timer = timer( ...
            'TimerFcn', @wd_callback, ...
            'Period', 5, ...
            'ExecutionMode', 'fixedRate' ...
        );
        start(this.wd_timer);
        
        function wd_callback(~,~)
            telapsed = toc(this.wd_timestamp);
            if telapsed > this.wd_keepalive_interval
                this.wd_timestamp = tic;
                this.keepalive();
            end
        end
    end

    function delete(this)
        stop(this.wd_timer);
        delete(this.wd_timer);
        m2c_close(this.handle);
    end

    function keepalive(this)
        disp('sending keepalive message');
        % send an abort message. This is the shortest possible 
        % command. It shouldn't hurt, because we only send keepalive
        % messages when the network connection is idle
        m2c_send_frame(this.handle,'%01#AB**');
        % if the abort message should cause problems, use the ones
        % below
        % m2c_send_frame(this.handle,'%01#RCP1X0000**');
        % m2c_recv_frame(this.handle);
    end
    
    function update_watchdog(this)
        this.wd_timestamp = tic;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function result = IEC_read(this, IEC_struct, varargin)
        this.update_watchdog();

        if nargin > 2
            IEC_struct.offset=varargin{1};
        end
        if nargin > 3
            IEC_struct.length=varargin{2};
        end
        if nargin > 4
          error('too many input arguments');
        end
        switch IEC_struct.location

          case 'I'
             result=this.IEC_readI(IEC_struct);
          case 'Q'
             result=this.IEC_readQ(IEC_struct);
          case 'M'
             result=this.IEC_readM(IEC_struct);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function result = IEC_write(this, IEC_struct, data, varargin)
        this.update_watchdog();

        if nargin > 3
            IEC_struct.offset=varargin{1};
        end
        if nargin > 4
            IEC_struct.length=varargin{2};
        else
            IEC_struct.length=length(data);
        end
        if nargin > 5
          error('too many input arguments');
        end
        switch IEC_struct.location

          % case 'I', makes no sense
          %   result=this.IEC_writeI(IEC_struct, data);
          case 'Q'
             result=this.IEC_writeQ(IEC_struct, data);
          case 'M'
             result=this.IEC_writeM(IEC_struct, data);
        end

    end

end

methods (Access = protected)
    
    function result = IEC_readI(this, IEC_struct)
        switch IEC_struct.type_id
          case 'X'
             result=this.IEC_readIX(IEC_struct);
          case 'W'
          case 'D'
       end
    end

    function result = IEC_readIX(this, IEC_struct)
        result = m2c_RC(this.handle, 1, 'X', IEC_struct.file, IEC_struct.element) == '1';
    end

    function result = IEC_readQ(this, IEC_struct)
        switch IEC_struct.type_id
          case 'X'
             result=this.readQX(IEC_struct);
          case 'W'
          case 'D'
       end
    end

    function result = IEC_readQX(this, IEC_struct)
        result = m2c_RC(this.handle, 1, 'Y', IEC_struct.file, IEC_struct.element) == '1';
    end

    function result = IEC_readM(this, IEC_struct)
        switch IEC_struct.type_id
          case 'X'
             result=this.IEC_readMX(IEC_struct);
          case 'W'
             result=this.IEC_readMW(IEC_struct);
          case 'D'
       end
    end

    function result = IEC_readMX(this, IEC_struct)
        result = m2c_RC(this.handle, 1, 'R', IEC_struct.element, IEC_struct.bit) == '1';
    end

    function result = IEC_readMW(this, IEC_struct)
        start_addr = IEC_struct.element+IEC_struct.offset;
        end_addr = start_addr + IEC_struct.length - 1;
        %result = m2c_RD(this.handle, 1, 'D', start_addr, end_addr);
        result = m2c_toint(m2c_RD(this.handle, 1, 'D', start_addr, end_addr));
    end

    function result = IEC_writeQ(this, IEC_struct, data)
        switch IEC_struct.type_id
          case 'X'
             result=this.writeQX(IEC_struct, data);
          case 'W'
          case 'D'
       end
    end

    function result = IEC_writeQX(this, IEC_struct, value)
        result = m2c_WC(this.handle, 1, 'Y', IEC_struct.file, IEC_struct.element, value);
    end

    function result = IEC_writeM(this, IEC_struct, data)
        switch IEC_struct.type_id
          case 'X'
             result=this.IEC_writeMX(IEC_struct, data);
          case 'W'
             result=this.IEC_writeMW(IEC_struct, data);
          case 'D'
       end
    end

    function result = IEC_writeMX(this, IEC_struct, value)
        result = m2c_WC(this.handle, 1, 'R', IEC_struct.element, IEC_struct.bit, value);
    end

    function result = IEC_writeMW(this, IEC_struct, data)
        start_addr = IEC_struct.element+IEC_struct.offset;
        %end_addr = start_addr + IEC_struct.length - 1;
        result = m2c_toint(m2c_WD(this.handle, 1, 'D', start_addr, data));
    end

end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function bResult = Read_X(this, addr, bit)
       this.update_watchdog();
       bResult = m2c_RC(this.handle, 1, 'X', addr, bit);
    end

    function bResult = Read_R(this, addr, bit)
       this.update_watchdog();
       bResult = m2c_RC(this.handle, 1, 'R', addr, bit);
    end

    function bResult = Read_Y(this, addr, bit)
       this.update_watchdog();
       bResult = m2c_RC(this.handle, 1, 'Y', addr, bit);
    end

    % !!!
    % function rResult = Read_Real(this, addr)
    % end

    % !!!
    % function Read_WX(this, addr)
    % end

    % !!!
    % function Read_WR(this, addr)
    % end

    % function Read_DWR(this, addr)
    % end

    % function Read_MultiWR(this, addr)
    % end

    % !!!
    % function Read_WY(this, addr)
    % end

    % !!!
    % function Read_Multi_WR_Bit(this, addr)
    % end

    % function ReadPLC_PV_DataTable(this,addr)
    % end

    % function ReadPLC_SV_DataTable(this,addr)
    % end

    % function Write_R(this, addr);
    % end

    % function Write_Real(this, addr);
    % end

end

end
