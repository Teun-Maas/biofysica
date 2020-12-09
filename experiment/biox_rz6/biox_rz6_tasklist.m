classdef biox_rz6_tasklist < handle
    properties (Access=protected)
        tq = [];        
        ip;
        nq = 0;
        debugflag = false;
    end
    
    properties (Access=public)

    end
    
    properties (Constant)
        desclen = 7;
        
        task_waitfortrigger = 0;
        task_sound_a = 1;
        task_sound_b = 2;
        task_sound_ab = 16;
        task_mux = 3;
        task_holdinput = 4; 
        task_soundmov = 5;
        task_daq = 6;
        task_dioout = 7;
        task_trgout = 8;
        task_reset = 9;
        task_ready = 10;
        task_multiconfiga = 11;
        task_multiconfigb = 12;
        task_att = 13;
        task_itd = 14;
        task_mix  = 15;
        
        soundtype_stop = 0;
        soundtype_tone = 1;
        soundtype_sweep = 2;
        soundtype_noise = 3;
        soundtype_ripple = 4;
        soundtype_wav = 5;
        soundtype_multitone = 6;
        soundtype_ba = 7;
    end

    methods
        function this = biox_rz6_tasklist 
        end
        
        function n = nr_of_tasks(this)
            n = this.nq;
        end       
               
        function tq_f32=get(this)
            tq_f32 = this.tq(1:this.nq,:);     
        end
        
        function add_task(this, delaytime, task, varargin)
           desc = this.parse_command(delaytime, task, varargin{:});
           this.append_task(desc);
        end

        function debug(this, val)
           if (val)
              this.debugflag = true;
           else
              this.debugflag = false;
           end
        end        
    end


    methods (Access=protected)
        function append_task(this, desc)
           this.nq = this.nq+1;
           this.tq(this.nq,:) = struct2array(desc);
        end

        function desc = parse_command(this, varargin)
            funcName='biox_rz6_tasklist.m/add_task';
            validCommands = { 'WaitForTrigger', 'SoundA','SoundB','Mux','HoldInput',...
               'SoundMov','Daq','DaqEx','SetDIO','TrigOut','Reset','Ready',...
               'MultiConfigA','MultiConfigB','Att','ITD','Mix','SoundAB'};

            validCommand = @(x) any(validatestring(x,validCommands));
            validDelayTime = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});

            p = biox_inputParser;
            p.FunctionName = funcName;
            p.addRequired('DelayTime', validDelayTime);
            p.addRequired('Command', validCommand);
            p.parse(varargin{1:2});

            % Expand partially matched commands 
            command = lower(validatestring(p.Results.Command,validCommands));
            switch command
                case 'waitfortrigger'
                    desc = this.parse_waitfortrigger(p,varargin{:});

                case 'sounda'
                    desc = this.parse_sound(this.task_sound_a,p,varargin{:});

                case 'soundb'
                    desc = this.parse_sound(this.task_sound_b,p,varargin{:});
                    
                case 'soundab'
                    desc = this.parse_sound(this.task_sound_ab,p,varargin{:});

                case 'mux'
                    desc = this.parse_mux(p,varargin{:});

                case 'holdinput'
                    desc = this.parse_holdinput(p,varargin{:});

                case 'soundmov'
                    desc = this.parse_soundmov(p,varargin{:});

                case 'daq'
                    desc = this.parse_daq(p,varargin{:});
                    
                case 'daqex'    
                    desc = this.parse_daqex(p,varargin{:});
                    
                case 'setdio'
                    desc = this.parse_setdio(p,varargin{:});

                case 'trigout'
                    desc = this.parse_trigout(p,varargin{:});

                case 'reset'
                    desc = this.parse_reset(p,varargin{:});

                case 'ready'
                    desc = this.parse_ready(p,varargin{:});

                case 'multiconfiga' %RL: naam veranderd
                    desc = this.parse_multiconfiga(p,varargin{:});

                case 'multiconfigb' %RL: naam veranderd
                    desc = this.parse_multiconfigb(p,varargin{:});

                case 'att'
                    desc = this.parse_att(p,varargin{:});

                case 'itd'
                    desc = this.parse_itd(p,varargin{:});
                    
                case 'mix'  %RL: toegevoegd   
                    desc = this.parse_mix(p,varargin{:});
                    
            otherwise
                error('unknown task: %s', p.Results.Command);
            end

            desc.DelayTime = p.Results.DelayTime;

            if this.debugflag
               disp(['---',mfilename,' debug output---']);
               disp('INPUT:');
               disp(p.Results);
               disp('');
               disp('RESULT:');
               disp(desc);
               disp('---end---');
            end

            if any(isnan(struct2array(desc)))
               error('missing arguments');
            end
        end

        function desc = newdesc(~)
            desc = struct( ...
              'TaskType', 0, ...
              'SoundType', 0, ...
              'DelayTime', 0, ...
              'Par1', 0, ...
              'Par2', 0, ...
              'Par3', 0, ...
              'Par4', 0 ...
            );
        end

        function desc = parse_waitfortrigger(this,p,varargin)
            desc = this.newdesc();
            desc.TaskType = this.task_waitfortrigger;

            validExtTrig = @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative','<',256}); %aangpast RL 8-->256            
                    
            expectedStringInputs = { 'ZBusA', 'ZBusB', 'External', 'Soft1', 'Soft2', 'Soft3'};  %RL: extra opties toegevoegd
            if ischar(varargin{3})
               isValid = @(x) any(validatestring(x,expectedStringInputs));               
            else
               isValid = @(x) validateattributes(x, {'numeric'},{'scalar','nonnegative','<',64}); %RL: 8 veranderd in 64.
            end            
            p.addRequired('Input', isValid);            
            
            p.addOptional('ExternalTrigger', 0, @(x) validExtTrig(x));
            p.parse(varargin{:});
            
            desc.Par1 = input2num(p.Results.Input);
            desc.Par2 = p.Results.ExternalTrigger;

           

            function n = input2num(x) 
               % expectedStringInputs = {  'ZBusA', 'ZBusB', 'External', 'Soft1', 'Soft2', 'Soft3' };
               if ischar(x)  
                  x = lower(x); 
                  if     x(5) == 'a'     %RL: extra opties toegevoegd
                     n = 1;
                  elseif x(5) == 'b'
                     n = 2;
                  elseif x(1) == 'e'
                     n = 4;
                  elseif x(5) == '1'
                     n = 8;                     
                  elseif x(5) == '2'
                     n = 16;                                          
                  elseif x(5) == '3'
                     n = 32;                                                               
                  end
               else % when x is a number
                  n = x;
               end
            end
        end

        function desc = parse_sound(this,tasktype,p,varargin)
            desc = this.newdesc();
            
            desc.TaskType = tasktype;

            expectedSounds = {'Stop','Tone','Sweep','Noise','Ripple','WAV','B=A','MultiTone'}; %RL 'MultiTone' toegevoegd

            validSound      = @(x) any(validatestring(x, expectedSounds));            
            validFreq       = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1e+006});
            validITD        = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
            validPosNum     = @(x) validateattributes(x,{'numeric'},{'scalar','positive'});            
            valid01         = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1});
            validStartPhase = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-180,'<=',180});

            p.addRequired('SoundType', validSound);
            p.parse(varargin{1:3});

            % Expand partially matched SoundType strings, and convert to lowercase
            SoundType=lower(validatestring(p.Results.SoundType, expectedSounds));
            
            function n = input2bool(x)   %RL: for use in 'noise' case
               if ischar(x)  
                  x = lower(x); 
                  if     x(1) == 's'   
                     n = 0;
                  elseif x(1) == 'c'
                     n = 1;
                  end
               end
            end   

            switch SoundType
            case 'stop'
               desc.SoundType = this.soundtype_stop;
               p.parse(varargin{:});

            case 'tone'
               desc.SoundType = this.soundtype_tone;
               p.addRequired('ToneFreq', validFreq);
               p.addOptional('ModFreq', 0, validFreq); %RL optional van gemaakt
               p.addOptional('ModBW', 0, validFreq);   %RL optional van gemaakt
               p.parse(varargin{:});
               desc.Par1 = p.Results.ToneFreq;
               desc.Par2 = p.Results.ModFreq;
               desc.Par3 = p.Results.ModBW;
               
            case 'multitone'
               desc.SoundType = this.soundtype_multitone;
               p.parse(varargin{:});               

            case 'sweep'
               desc.SoundType = this.soundtype_sweep;
               p.addRequired('StartFreq', validFreq);
               p.addRequired('NrOctaves', validPosNum);
               p.addRequired('Period', validPosNum);
               p.parse(varargin{:});
               desc.Par1 = p.Results.StartFreq;
               desc.Par2 = p.Results.NrOctaves;
               desc.Par3 = p.Results.Period; % Period in msec

            case 'noise'               
               desc.SoundType = this.soundtype_noise;
               p.addRequired('HpFreq1', validFreq);
               p.addRequired('LpFreq1', validFreq);
              
               p.parse(varargin{1:5});
               
               HpFreq1 = p.Results.HpFreq1;
               LpFreq1 = p.Results.LpFreq1;  
                             
               p.addOptional('HpFreq2', HpFreq1, validFreq);
               p.addOptional('LpFreq2', LpFreq1, validFreq);
               
               p.parse(varargin{:});
               
               HpFreq2 = p.Results.HpFreq2;
               LpFreq2 = p.Results.LpFreq2;                                 
                   
               if HpFreq1 >= LpFreq1
                 error('LpFreq must be greater than HpFreq');  
               end 
               
               if (HpFreq2 >= LpFreq2) && (HpFreq ~= 0)
                 error('LpFreq must be greater than HpFreq');  
               end 
               
               desc.Par1 = HpFreq1;
               desc.Par2 = LpFreq1;
               desc.Par3 = HpFreq2;
               desc.Par4 = LpFreq2;

            case 'ripple'
               desc.SoundType = this.soundtype_ripple;
               p.addRequired('StartFreq', validFreq);
               p.addRequired('ModInTime', validPosNum);
               p.addRequired('ModInFreq', validPosNum);
               p.addRequired('ModDepth' , valid01);
               p.parse(varargin{:});
               desc.Par1 = p.Results.StartFreq;
               desc.Par2 = p.Results.ModInTime; % Hz
               desc.Par3 = p.Results.ModInFreq; % Phase/Octave
               desc.Par4 = p.Results.ModDepth; % Modulation depth

            case 'wav'
               desc.SoundType = this.soundtype_wav;
               %RL mischien 'Reset' als input ipv 1?
               expectedPar1 = { 'Reset', 'Continue', 'Loop'};
               p.addOptional('Reset', 'Reset', @(x) any(validatestring(x, expectedPar1))); %RL optional ipv Required
               p.parse(varargin{:});
               switch lower(p.Results.Reset)
               case 'continue'
                  desc.Par1 = 0;  
               case 'reset'
                  desc.Par1 = 1;
               case 'loop'   
                  desc.Par1 = 2; 
               otherwise
                 error('invalid parameter for WAV-sound, this is a bug');                      
               end  

            case 'b=a'
               assert(desc.TaskType == this.task_sound_b,...
                  'SoundType ''b=a'' is only valid for taskType ''SoundB''');
               desc.SoundType = this.soundtype_ba;
               expectedMovTypes = { 'Fixed','Linear','Sine' };
               p.addRequired('movType', @(x) any(validatestring(x, expectedMovTypes)));
               p.parse(varargin{1:4});
               movType=lower(validatestring(p.Results.movType,expectedMovTypes));
               switch lower(movType)
               case 'fixed'
                  p.parse(varargin{:});
                  desc.Par1 = 1; %aangepast RL
               case 'sine'
                  p.addRequired('Period', validPosNum);
                  p.addRequired('Phase', validStartPhase);
                  p.parse(varargin{:});
                  desc.Par1 = 2; %aangepast RL
                  desc.Par2 = p.Results.Period; % msec
                  desc.Par3 = p.Results.Phase; 
               case 'linear'
                  p.addRequired('Period', validPosNum);
                  p.addRequired('Phase', validStartPhase);
                  p.parse(varargin{:});
                  desc.Par1 = 3;  %aangepast RL
                  desc.Par2 = p.Results.Period; % msec
                  desc.Par3 = p.Results.Phase; 
               otherwise
                  error('invalid mov type, this is a bug');
               end

            otherwise
               error('unknown sound type, this is a bug');
            end

        end

        function desc = parse_mux(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_mux;

           validChan = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',15});
           validDevice = @(x) validateattributes(x, {'numeric'}, {'scalar','nonnegative','<=',3});

           expectedActions = { 'Set','Reset' };
           p.addRequired('Device', @(x) validDevice(x));
           p.addRequired('Action', @(x) any(validatestring(x, expectedActions)));
           p.addOptional('Channel', [], @(x) validChan(x));

           p.parse(varargin{:});

           if strcmpi(p.Results.Action,'Set')
               if isempty(p.Results.Channel)
                   error('''Set'' requires a the value of ''Channel'' to be specified.');
               else
                   Channel = p.Results.Channel;
               end
               ASet=1; AReset=0;
           elseif strcmpi(p.Results.Action,'Reset')
               Channel = 0;
               ASet=0; AReset=1;
           else
               error('Unexpected or missing ''Action'', this is a bug');
           end   
           desc.Par1 = Channel + 16*p.Results.Device + 64*ASet + 128*AReset;
        end

        function desc = parse_holdinput(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_holdinput;

           validInputMask = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative', '<',256});

           p.addRequired('InputMask', validInputMask);
           p.parse(varargin{:});
           desc.Par1 = p.Results.InputMask;
        end

        function desc = parse_soundmov(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_soundmov;

           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));
           validNumSpeakers = @(x) validateattributes(x, ...
              {'numeric'}, {'odd','scalar','>=',3,'<=',21});
           validPeriod = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative'});
           validStartPhase = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>=',-180, '<=',180});

           p.addRequired('StartStop', validStartStop);
           p.parse(varargin{1:3});
           switch lower(p.Results.StartStop)
           case 'stop'
               p.parse(varargin{:});
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
               p.addRequired('NumSpeakers', validNumSpeakers);
               p.addRequired('Period', validPeriod);
               p.addRequired('StartPhase', validStartPhase);
               p.parse(varargin{:});
               desc.Par2 = p.Results.NumSpeakers;
               desc.Par3 = p.Results.Period;
               desc.Par4 = p.Results.StartPhase;
           otherwise
               error('unexcpected error, this is a bug');
           end
        end

        function desc = parse_daq(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_daq;
           validChannelSelection = @(x) validateattributes(x, ...
               {'numeric'}, {'vector','nonnegative'});
           validDivisor = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','positive'});
           validStartStop = @(x) any(validatestring(x,{'Start','Stop'}));
         
           p.addRequired('StartStop', validStartStop);
           p.addRequired('ChannelSelection', validChannelSelection);           
           p.addOptional('Divisor', 1, validDivisor);

           p.parse(varargin{:});
           switch lower(p.Results.StartStop)
           case 'stop'
               desc.Par1 = 0;
           case 'start'
               desc.Par1 = 1;
           otherwise
               error('unexcpected error, this is a bug');
           end
          
           ChannelList = p.Results.ChannelSelection; % RL
           ChannelListInt = sum(2.^ChannelList); %RL: convert to integer 
           
           desc.Par2 = ChannelListInt;
           desc.Par3 = round(p.Results.Divisor); %must be an integer
        end               

        function desc = parse_setdio(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_dioout;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});

           p.addRequired('OutputByte', validByte);
           p.parse(varargin{:});
           desc.Par1 = p.Results.OutputByte;
        end

        function desc = parse_trigout(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_trgout;

           validByte = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','nonnegative','<',256});
           validDTInterval = @(x) validateattributes(x, ...
              {'numeric'}, {'scalar','>=',0.000040}); %RL: changed to >= 0.000040 = two clock tics

           p.addRequired('OutputByte', validByte);                     
           p.addOptional('DoubleDelay',0, validDTInterval);           
           p.parse(varargin{:});           
           desc.Par1 = p.Results.OutputByte;          
           desc.Par2 = p.Results.DoubleDelay;
        end

        function desc = parse_reset(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_reset;
           p.parse(varargin{:});
        end

        function desc = parse_ready(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_ready;
           p.parse(varargin{:});
        end

        function desc = parse_multiconfig(this,tasktype,p,varargin)
           validIndex = @(x) validateattributes(x,{'numeric'},{'scalar','positive','<=',4});
           validFrequency = @(x) validateattributes(x,{'numeric'},{'scalar','nonnegative'});
           validPhase = @(x) validateattributes(x,{'numeric'},{'scalar','>=',-180,'<=',180});
           validAmplitude = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1000}); %RL bereik veranderd 0...1000%

           desc = this.newdesc();
           desc.TaskType = tasktype;
           p.addRequired('Index', validIndex);
           p.addRequired('Frequency', validFrequency);
           p.addRequired('Phase', validPhase);
           p.addRequired('Amplitude', validAmplitude);
           p.parse(varargin{:});
           desc.Par1 = p.Results.Index;
           desc.Par2 = p.Results.Frequency;
           desc.Par3 = p.Results.Phase;
           desc.Par4 = p.Results.Amplitude;
        end

        function desc = parse_multiconfiga(this,p,varargin)
            desc = this.parse_multiconfig(this.task_multiconfiga, p, varargin{:});
        end

        function desc = parse_multiconfigb(this,p,varargin)
            desc = this.parse_multiconfig(this.task_multiconfigb, p, varargin{:});
        end

        function desc = parse_att(this,p,varargin)
           validAttenuation = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',80});
           validScaleFactor = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1000});
           desc = this.newdesc();
           desc.TaskType = this.task_att;
           p.addRequired('Attenuation1', validAttenuation);
           p.addRequired('Attenuation2', validAttenuation);
           p.addOptional('ScaleFactor1', 1, validScaleFactor);
           p.addOptional('ScaleFactor2', 1, validScaleFactor);
           p.parse(varargin{:});
           desc.Par1 = p.Results.Attenuation1;
           desc.Par2 = p.Results.Attenuation2;
           desc.Par3 = p.Results.ScaleFactor1;
           desc.Par4 = p.Results.ScaleFactor2;
        end 

        function desc = parse_itd(this,p,varargin)
           validITD = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',0.01});           
           desc = this.newdesc();
           desc.TaskType = this.task_itd;
           p.addRequired('ITD1', validITD);
           p.addRequired('ITD2', validITD);
           p.parse(varargin{:});
           desc.Par1 = p.Results.ITD1;
           desc.Par2 = p.Results.ITD2;
        end         

         %RL functie parse_mix toegevoegd
        function desc = parse_mix(this,p,varargin)
           desc = this.newdesc();
           desc.TaskType = this.task_mix;

           validType =  @(x) any(validatestring(x,{'Stop','BtoA','AtoB','Mixed'}));
           
           validFactor = @(x) validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1});
           
           p.addRequired('Input', validType);                     
           p.addOptional('Factor',1,validFactor); 
           
           p.parse(varargin{:});           
           desc.Par1 = input2num(p.Results.Input);          
           desc.Par2 = p.Results.Factor;
           
           function n = input2num(x) 
           % expectedStringInputs = {'Stop','BtoA','AtoB','MixToBoth'};
               if ischar(x)  
                  x = lower(x); 
                  if     x(1) == 's'     %RL: extra opties toegevoegd
                     n = 0;
                  elseif x(1) == 'b'
                     n = 1;
                  elseif x(1) == 'a'
                     n = 2;
                  elseif x(1) == 'm'
                     n = 3;                                   
                  end
               else
                  n = 0;
               end
           end
        end    

    end

    methods (Static)

        function desc=define_task(task_type, sound_type, delay_time, par1, par2, par3, par4)
           desc = [ task_type sound_type delay_time par1 par2 par3 par4 ];
        end
        
        
    end
end
