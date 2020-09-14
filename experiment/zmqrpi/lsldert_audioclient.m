classdef lsldert_audioclient < lsldert_pubclient
    properties
        maxbuf = 10;
    end
    
    methods
        function this = lsldert_audioclient(varargin)
            % LSLDERT_AUDIOCLIENT class constructor
            %
            % obj = lsldert_audioclient('lslsdert-host.local');
            % creates a connection to the lsldert server listening
            % on port 5555 at remote host 'lslsdert-host.local'.
            % the port number can optionally be specified as a 2nd argument
            % obj = lsldert_client('lsldert-host.local', 5555);         
            this@lsldert_pubclient(varargin{:});
        end
        
        function load(this,bufferno, varargin)
            if (bufferno < 0) || (bufferno > this.maxbuf)
                ME = MException('lsldert_audioclient.load',...
                    'bufferno must be >=0 and <%d',...
                    this.maxbuf);
                throw(ME);
            end
            [pcm_data,Fs]=audioread(varargin{:});
            if ~(Fs==44100 || Fs==48000)
                Fs2=48000;
                fprintf('lsldert_audioclient: resampling data from %d to %d S/s\n',Fs,Fs2);
                pcm_data2=resample(pcm_data,Fs2,Fs);
                pcm_data=pcm_data2;
                Fs=Fs2;
            end
            nsamp=max(size(pcm_data));
            nchan=min(size(pcm_data));
            pcm_header=uint32([Fs, nchan, nsamp, 0, 0, 0]);
            cmd=sprintf('AF %d',bufferno);
            this.send(cmd, pcm_header, single(pcm_data));
        end
        
        function play(this,bufferno)
            if (bufferno < 0) || (bufferno > this.maxbuf)
                ME = MException('lsldert_audioclient.play',...
                    'bufferno must be >=0 and <%d',...
                    this.maxbuf);
                throw(ME);
            end
            cmd=sprintf('AP %d',bufferno);
            this.send(cmd);
        end
        
        function stop(this,bufferno)
            if (bufferno < 0) || (bufferno > this.maxbuf)
                ME = MException('lsldert_audioclient.stop',...
                    'bufferno must be >=0 and <%d',...
                    this.maxbuf);
                throw(ME);
            end
            cmd=sprintf('AS %d',bufferno);
            this.send(cmd);
        end
    end
end
