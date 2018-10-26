classdef ledcontroller_pi < handle
    % LEDCONTROLLER_PI is a class to program the led array on Raspberry
    % PI/Arduino modules. 
    %
    % See also LEDCONTROLLER_PI, DELETE, TRIGGER, TRIGGER_ENABLE, TRIGGER_COUNT, TRIGGERED,
    % WRITE, WAIT
    %
    % We need jeromq here! Add the full path to jeromq.jar to you javaclasspath.txt

    
    % Written by Günter Windau
    % 2014-04-22 version 1.0
    % 2014-04-22 version 1.1 optimized some parts manipulating individual bits
    % 2014-06-02 version 1.2 put led stimuli in a separate class, all stimuli must be triggered
    % 2014-06-02 version 1.3 major code cleanup, adaptation to improved
    % PLC program
    % 2017-10-19 version 2.0 adaptation for PI/Arduino hardware. Roughly
    % the same API as for the PLC version
    % 2018-01-16 GW: fixed intensity problem in write() function

    
    properties (Access=protected)
        hostnames = {''};
        context;
        sockets;
        nbuf;  % nr of buffers in PLC for playing LED stimuli.
        nbuf_written = 0; % nr of buffers written in ledpattern/write
        default_timeout = 30;
    end
    
    methods
        function this = ledcontroller_pi(varargin)
            % LEDCONTROLLER_PI(hostname [, hostname2 [,...])
            % Setup a connection to the ArduinoPI led driver at internet addresses or
            % hostnames specified
            import org.zeromq.*

            if isempty(varargin)
                error('no hosts specified');
            end
            this.hostnames = varargin;
                        
            this.context = zmq.Ctx();
            this.sockets = cell(size(this.hostnames));
            i=1;
            for h=this.hostnames
                this.sockets{i} = this.context.createSocket(ZMQ.REQ);
                this.sockets{i}.connect(strcat('tcp://', h, ':5555'));
            %   TODO 
            %   this.sockets{i}.setSocketOpt(ZMQ.RCVTIMEO,1000); FIXME
                i=1+1;
            end
            %this.print_version;
            this.nbuf=256;
            this.trigger_enable(1);
            this.clear_all;
        end
        
        function print_version(this)
            % PRINT_VERSION
            % print version information for the connected Arduino modules
            i=1;
            versions=this.send_command('RV');
            for h=this.hostnames
                hn=string(h);
                v=string(versions{i});
                fprintf('%s: %s\n', hn, v);
                i=i+1;
            end
        end    
            
        function delete(this)
            % DELETE - class destructor
            % Mop up and delete the LEDCONTROLLER_PI thisect
            cellfun(@(x) x.close, this.sockets);
            clear this.context;
        end

    end
    methods(Access = protected)
        function cells=split_args(this, args, num_cells, args_per_cell)
           elements_needed=num_cells*args_per_cell;
           ra=reshape(args,1,[]);  % make sure we have a row vector here
           
           if length(ra) < elements_needed
               % pad with zeros to have multiples of 8 args
               ra=[ra zeros(1,elements_needed-length(ra))];
           else
               % truncate 
               ra=ra(1:elements_needed);
           end
           % reshape ra 
           ra=transpose(reshape(ra, [args_per_cell num_cells]));
           cells=mat2cell(ra, ones(1,num_cells))';
        end
        
        function result=send_command(this, command, args)
            % SEND_COMMAND - command routing to clusters of AduinoPI led
            % drivers
            single_wr_commands = {'ER','WE','WT','WC','WI'}; % Write commands with one argument routed to all units
            single_rd_commands = {'T','RE','RT','RC'};  % Read commands expection one result from each unit
            multi_wr_commands = {'WR','WG'};  % Write commands with array arguments to be split across units
            multi_rd_commands = {'RV','RR','RG','RI','RJ'};  % Read commands returning array results to be combined from all units
    
         %   command
         %   if nargin>2 
         %       args
         %   end
            if any(ismember(single_wr_commands, command))
               cellfun(@(s) this.send_command_to(s, command, args), this.sockets);
               r=cellfun(@(s) this.receive_from(s), this.sockets, 'UniformOutput', 0);
               result=cell2mat(r);
            
            elseif any(ismember(single_rd_commands, command))
               cellfun(@(s) this.send_command_to(s, command), this.sockets);
               r=cellfun(@(s) this.receive_from(s), this.sockets, 'UniformOutput', 0);
               result=cell2mat(r);

            elseif any(ismember(multi_wr_commands, command))
               num_sockets=numel(this.sockets);
               ncol=this.nbuf;
               ca=this.split_args(args, num_sockets, ncol);
               socks=this.sockets;
               cellfun(@(s,a) this.send_command_to(s, command, a), socks, ca);
               result=cellfun(@(s) this.receive_from(s), this.sockets, 'UniformOutput', 0);
            
            elseif any(ismember(multi_rd_commands, command))
               cellfun(@(s) this.send_command_to(s, command), this.sockets);
               result=cellfun(@(s) this.receive_from(s), this.sockets, 'UniformOutput', 0);
            end
        end

        function send_command_to(~, socket, command, args)
            import org.zeromq.*;

            if nargin <= 3
                cmd = command;
            else
                cmd=strcat(command, sprintf(' %d', args), char(0));
            end
            
            %fprintf('send_command sent: %s\n', cmd);
            message = zmq.Msg(1024);
            message.put(unicode2native(cmd));

            socket.send(message,0);
        end

        function result = receive_from(this, socket)
            nonblock=1;
            message=[];
            t0=tic;
            while 1  % Need to find out how to set ZMQ_RCVTIMEO here...
                message=socket.recv(nonblock);
                if ~isempty(message)
                    break;
                end   
                pause(0.01);
                if toc(t0)>20
                    error('recv timeout');
                end
            end
            r=native2unicode(message.data)';
             %fprintf('receive_from recv: %s\n', r);

            if (all(r(1:2)=='ER') || all(r(1:2)=='RX'))
                  % get verbose error message first, then call error()
                  r(1:2)='RX';
                  this.send_command_to(socket,r);
                  message=socket.recv(0);
                  r=native2unicode(message.data)';
                  r=r(4:end);
                  error('send_command: %s', r);   
            elseif r(1:2)=='RV'
                result=r(4:end);
            elseif r(1)=='R'
               result=str2num(r(4:end));
            elseif r(1:2)=='OK'
               result=[];
            else
               error('send_command received unexpected answer: %s', r);
            end
        end
     
    end
    
    methods (Access = public)
        function trigger_enable(this, value)
            % TRIGGER_ENABLE(value)
            % Set or clear the trigger enable bit in the control register.
            % If the bit is set, the PLC program will flush the buffer to
            % the outputs when a rising edge trigger is detected on the TTL
            % trigger input
            this.send_command('WE', value);
        end
        
        function trigger_wait(this, timeout)
            % TRIGGER_WAIT
            if nargin < 2
                timeout = this.default_timeout;
            end
            t0=tic;
            while ~this.triggered
                if toc(t0) > timeout
                    error('trigger_wait timed out');
                end
                % not cpu friendly, but for now loop around
            end
        end
       
        function result = trigger_count(this)
            % TRIGGER_COUNT
            % Return the number of trigger since (?)
            result = this.send_command('RC');
        end
        
        function trigger(this)
            % TRIGGER
            % Trigger by software
            %tic;
            this.send_command('T');
            %deltat=toc
        end
        
        function result = triggered(this)
            % TRIGGERED
            % Return true if a rising edge was detected on the external
            % trigger input.
            trigd = this.send_command('RT');
            if trigd
                this.send_command('WT', 0);
            end
            result = trigd;
        end
        
        function gate_enable(this, value)
            % GATE_ENABLE(value)
            % set or clear the gate enable bit in the control register.
            % If the bit is set, the PLC program will flush the buffer to
            % the outputs immediately when the TTL gate input is at TTL
            % high level.
            % not implemented in raspi/arduino hardware
            % set/reset bit nr. 7
            warning('ledcontroller_pi.gate_enable() is not implemented');
        end
        
             
        function clear_all(this)
            % CLEAR_ALL
            % Turn all leds off. This function is optimized, and changes
            % are written to the PCL immediately, but still depending on
            % trigger and gate signals.

            this.send_command('WR', [0, 0, 0, 0, 0, 0, 0, 0]);
            this.send_command('WG', [0, 0, 0, 0, 0, 0, 0, 0]);
        end
        
        function write(this, stimulus)
            % WRITE
            % write LED arrays to the PLC. 
            % The outputs change after each trigger.
            this.trigger_enable(0);
            [m,n] = size(stimulus);
            num_sockets=numel(this.sockets);
            
            if m*n>this.nbuf
               error('ledcontroller_pi/write: too many buffers: %d>%d', m*n, this.nbuf); 
            end
            
            ibuf = 1;
            all_red = zeros(this.nbuf, num_sockets, 'uint16');
            all_grn = zeros(this.nbuf, num_sockets, 'uint16');
            all_i = zeros(this.nbuf, 1, 'uint16');
            for i = 1:m
                for j = 1:n
                    [leds_red, leds_grn, ir, ig] = stimulus(i,j).get_leds;
                    
                    %Handle intensities
                    %because we have only one intesity control on arduino
                    %let's clean up de red, green and intensity settings
                    %in order to get the best result
                    if any(leds_red) 
                        if (ir==0)
                           leds_red(:)=0;
                        end
                    else
                        ir=0;
                    end
                    if any(leds_grn) 
                        if (ig==0)
                           leds_grn(:)=0;
                        end
                    else
                        ig=0;
                    end
                    if (ir > 0) && (ig > 0) && (ir ~= ig)
                        warning('(ir,ig)=(%d,%d),Arduino LED driver cannot set different intensities for red and green',ir,ig);
                    end
                    intens=max(ir,ig);
                    if intens>50
                        warning('LED Intensity must be <= 50')
                        intens=50;
                    elseif intens < 0
                        warning('LED Intensity must be >= 0')
                        intens=0;
                    end
                    all_i(ibuf) = uint16(255.0*intens/50);

                    %Handle LED bits
                    %pack all boolean values into 16 bit words
                    T = 2.^(0:15)';  % bit transformation matrix
                 %BUG   red = reshape(leds_red, 16, this.nbuf)';
                 %BUG   grn = reshape(leds_grn, 16, this.nbuf)';
                    red = reshape(leds_red, 16, [])';
                    grn = reshape(leds_grn, 16, [])';
                    red_bits = uint16(red*T);
                    grn_bits = uint16(grn*T);
                    
                    all_red(ibuf,:) = red_bits(1:num_sockets);
                    all_grn(ibuf,:) = grn_bits(1:num_sockets);
                                  
                    ibuf=ibuf+1;
                end
            end
            all_red=reshape(all_red, [1 numel(all_red)]);
            all_grn=reshape(all_grn, [1 numel(all_grn)]);
            this.send_command('WR', all_red);
            this.send_command('WG', all_grn);
            this.send_command('WI', all_i);
            this.send_command('WC', 0);  % clear trigger count
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
    end
    
 end
