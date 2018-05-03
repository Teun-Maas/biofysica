classdef ledcontroller < handle
    % LEDCONTROLLER is a class to program the led array on Panasonic PLCs using modbus
    % methods
    %
    % See also LEDCONTROLLER, DELETE, TRIGGER_ENABLE, TRIGGER_COUNT, TRIGGERED, GATE_ENABLE,
    % WRITE, WAIT, DUMP_REGISTERS
    
    % Written by Günter Windau
    % 2014-04-22 version 1.0
    % 2014-04-22 version 1.1 optimized some parts manipulating individual bits
    % 2014-06-02 version 1.2 put led stimuli in a separate class, all stimuli must be triggered
    % 2014-06-02 version 1.3 major code cleanup, adaptation to improved
    % PLC program

    properties (Access=protected)
        hostname  = 'fp-web2';
        baseaddr  = 55000;
        ctlbits   = 55000;
        trigcnt   = 55001;
        stimttl  = 55002; % LED Time to Live, in PLC cycles, 0=Inf (default)
        basered   = 55100;
        basegrn   = 55108;
        pwm0dc    = 55022; % red, duty cycle, 32 bits, range 0..100
        pwm1dc    = 55026; % green
        pwm2dc    = 55030; % red
        pwm3dc    = 55034; % green
        baseint   = 55228; % intensity per buffer: these are four groups of eight uint32's (8 buffers for each PWM)
        basenbuf  = 55292;
        ctlbits_trigger_enable = 9;
        ctlbits_gate_enable = 8;
        ctlbits_triggered = 7;
        
        hmodbus;
        cache_ctlbits;
        nbuf;  % nr of buffers in PLC for playing LED stimuli.
        nbuf_written = 0; % nr of buffers written in ledcontroller/write
    end
    
    methods
        function obj = ledcontroller(hostname)
            % LEDCONTROLLER(hostname)
            % Setup a connection to the PLC at internet address or
            % hostname specified
            if nargin > 0
                obj.hostname = hostname;
            end
            ip = obj.resolveip(obj.hostname);
            obj.hmodbus = modbus_new_tcp(ip, 502);
            modbus_connect(obj.hmodbus);
            
            obj.nbuf=modbus_read_registers(obj.hmodbus, obj.basenbuf, 1);
            if obj.nbuf <= 0
               error('ledcontroller/ledcontroller: according to the PLC the number of LED buffers is %d.', obj.nbuf);
               % obj.nbuf = 8;
            end
            % Enable the PWM (bits 3-0), trigger is enabled, and gate is disabled
            obj.cache_ctlbits = uint16(15);
            modbus_write_registers(obj.hmodbus, obj.ctlbits, obj.cache_ctlbits);
            
            obj.trigger_enable(1);

            obj.clear_all;
        end
        
        function delete(obj)
            % DELETE - class destructor
            % Mop up and delete the LEDCONTROLLER object
            obj.intensity('r',50);
            obj.intensity('g',50);
            modbus_close(obj.hmodbus);
            modbus_free(obj.hmodbus);
        end
        
        function set_stimttl(this,count)
            % SET_STIMTTL(count)
            % Set the time-to-live for a stimulus in PLC cycles.
            % This allows for very short led flashes.
            % By default count=0, LED's stay on until changed by the
            % next stimulus. 
            modbus_write_registers(this.hmodbus, this.stimttl, count);
        end

        function trigger_enable(this, value)
            % TRIGGER_ENABLE(value)
            % Set or clear the trigger enable bit in the control register.
            % If the bit is set, the PLC program will flush the buffer to
            % the outputs when a rising edge trigger is detected on the TTL
            % trigger input
            this.cache_ctlbits = bitset(this.cache_ctlbits, this.ctlbits_trigger_enable, value);
            modbus_write_registers(this.hmodbus, this.ctlbits, this.cache_ctlbits);
        end
        
        function trigger_wait(this)
            % TRIGGER_WAIT
            while ~this.triggered
                % not cpu friendly, but for now loop around
            end
        end
       
        function result = trigger_count(this)
            % TRIGGER_COUNT
            % Return the number of trigger since (?)
            result = modbus_read_registers(this.hmodbus, this.trigcnt, 1);
        end
        
        function result = triggered(this)
            % TRIGGERED
            % Return true if a rising edge was detected on the external
            % trigger input.
            reg = modbus_read_registers(this.hmodbus, this.ctlbits, 1);
            trigbit = bitget(reg, this.ctlbits_triggered);
            if trigbit > 0
                reg = uint16(bitset(reg, this.ctlbits_triggered, 0)); % clear triggered bit
                modbus_write_registers(this.hmodbus, this.ctlbits, reg);
                result = true;
            else
                result = false;
            end
        end
        
        function gate_enable(this, value)
            % GATE_ENABLE(value)
            % set or clear the gate enable bit in the control register.
            % If the bit is set, the PLC program will flush the buffer to
            % the outputs immediately when the TTL gate input is at TTL
            % high level.
            % set/reset bit nr. 7
            this.cache_ctlbits = uint16(bitset(this.cache_ctlbits, this.ctlbits_gate_enable, value));
            modbus_write_registers(this.hmodbus, this.ctlbits, this.cache_ctlbits);
        end
        
             
        function clear_all(this)
            % CLEAR_ALL
            % Turn all leds off. This function is optimized, and changes
            % are written to the PCL immediately, but still depending on
            % trigger and gate signals.
            bits = zeros(8, 1, 'uint16');
            modbus_write_registers(this.hmodbus, this.basered, bits);
            modbus_write_registers(this.hmodbus, this.basegrn, bits);
        end
        
        function write(this, stimulus)
            % WRITE
            % write LED arrays to the PLC. 
            % The outputs change after each trigger.
            this.trigger_enable(0);
            pwm_dc = zeros(this.nbuf*2, 4, 'uint16'); % PWM Duty Cycle map
            [m,n] = size(stimulus);
            
            if m*n>this.nbuf
               error('ledcontroller/write: too many buffers: %d>%d', m*n, this.nbuf); 
            end
            
            ibuf = 1;
            addr = this.basered;
            for i = 1:m
                for j = 1:n
                    [leds_red, leds_grn, ir, ig] = stimulus(i,j).get_leds;
                    
                    % pack all boolean values into 16 bit words
                    T = 2.^[0:15]';  % bit transformation matrix
                    red = reshape(leds_red, 16, 8)';
                    grn = reshape(leds_grn, 16, 8)';
            
                    red_bits = uint16(red*T);
                    grn_bits = uint16(grn*T);
            
                    %modbus_write_registers(this.hmodbus, this.basered, red_bits);
                    %modbus_write_registers(this.hmodbus, this.basegrn, grn_bits);
                    % the above two lines optimized into one modbus write operation
                    all_bits = [ red_bits; grn_bits ];
                    modbus_write_registers(this.hmodbus, addr, all_bits);
                    addr = addr + 16;

                    % Set the intensity of the red or green leds
                    % if color starts with 'r' the intensity of the red leds is set
                    % by manpulating the PWM output duty cycle. Otherwise the green
                    % intensity is set.
                    % The value should be in the range 1..50
                    % pack all intensity values into arrays
                    pwm_dc(2*(ibuf-1)+1, [1 3]) = uint16(ir);
                    pwm_dc(2*(ibuf-1)+1, [2 4]) = uint16(ig);
                    
                    ibuf=ibuf+1;
                end
            end
            % write duty cylces to PLC
            %fprintf('pwm_dc=\n'); disp(pwm_dc);
            pwm_dc=reshape(pwm_dc,2*this.nbuf*4,1);
            %fprintf('pwm_dc=\n'); disp(pwm_dc);
            modbus_write_registers(this.hmodbus, this.baseint, pwm_dc);
            modbus_write_registers(this.hmodbus, this.trigcnt, uint16(0));
            this.trigger_enable(1);
            this.nbuf_written=n*m;
        end
        
        function done=wait(this, timeout)
            % WAIT
            % WAIT(timeout)
            % wait for the buffers to complete on the PLC
            % but no longer than timeout seconds.
            % WAIT returns 1 if the operation is complete, or 0 if the timeout period was exceeded.
            %
            % If timeout is not specified WAIT will wait indefinitely.
            tstart=tic;

            if nargin < 2
                timeout=Inf;
            end
            
            while timeout > toc(tstart)
                count=this.trigger_count;
                %fprintf('trigger count = %d\n', count);
                if count >= this.nbuf_written
                   done=1;
                   return
                end
                pause(0.01); % be friendly to the CPU
            end            
            done=0;
        end
        
        function dump_registers(this)
            % DUMP_REGISTERS
            % Write a dump of all PLC registers to stdout. This function
            % may be useful for debugging.
            nregs = 35;
            regs = modbus_read_registers(this.hmodbus, this.baseaddr, nregs);
            addr = this.baseaddr;
            disp('ledcontroller/dump_registers: dump of PLC registers follows:');
            for i=1:nregs
                disp([' ' num2str(addr) ' (0x' dec2hex(addr) '): ' dec2bin(regs(i), 16)]);
                addr = addr+1;
            end
            disp('end.');
        end
    end
    
    methods (Access=protected)
        function intensity(this, color, value)
            % INTENSITY(color, value)
            % Set the intensity of the red or green leds
            % if color starts with 'r' the intensity of the red leds is set
            % by manpulating the PWM output duty cycle. Otherwise the green
            % intensity is set.
            % The value should be in the range 1..50
            if color(1) == 'r'
                modbus_write_registers(this.hmodbus, this.pwm0dc, uint16(value));
                modbus_write_registers(this.hmodbus, this.pwm2dc, uint16(value));
            else
                modbus_write_registers(this.hmodbus, this.pwm1dc, uint16(value));
                modbus_write_registers(this.hmodbus, this.pwm3dc, uint16(value));
            end
        end
        
   
        
    end
    methods(Static)
        function ipaddress=resolveip(input)
            try
                address = java.net.InetAddress.getByName(input);
                ipaddress = char(address.getHostAddress);
                %  hostname = char(address.getHostName);
            catch
                error('ledcontroller/resolveip: unknown host %s.', input);
            end
        end
        
    end
end
