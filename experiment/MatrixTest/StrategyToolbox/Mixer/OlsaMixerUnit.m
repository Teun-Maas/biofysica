% OlsaMixerUnit < ProcUnit
% Mixes two input audio signals (1:speech, 2:noise). Input 2 can be empty.  
% Target output levels must be specified. The speech signal is scaled based on a 
% user-specified nominal input rms if provided (i.e. non-empty). The target speech
% level is expressed in dB SPL rms. The noise level type needs to be 
% specified (dB rms/peak, or relative to the noise input). sensIn defines 
% the peak SPL equivalent to 0dB FS. wrap controls the behaviour if input 
% signals have different lengths warp-around / zero-pad noise to match 
% duration of the speech signal). For wrap=1, durFade controls the duration
% of a cosine-shaped cyclic cross-fade. noiseSeek determines the starting 
% position within the noise file (set to -1 for random position). 
%
% OlsaMixerUnit Properties:
%   nomSpeechRms - nominal speech input RMS on linear scale (can be empty)
%   sensIn    - input sensitivity (dB SPL equivalent to sinusoid at digital FS)
%   lvlType   - 'rms' / 'peak' / 'rel'; 
%   lvlDb     - 2-element vector of levels  in dB
%   delays    - vector of onset delays for each input [s] 
%   wrap      - repeat shorter inputs to match duration of the longest
%               input? [1/0]
%   durFade   - duration of cross-fade when wrapping signals [s]
%   noiseSeek - initial seek position within noise [s]; -1 for random
%   clip      - read-only indicator: did clipping occur during the most recent
%               call of run()? [0/1]
%
% AudioMixerUnit Methods:
%    OlsaMixerUnit - constructor
%    run - execute processing
%
% See also: olsaMixerFunc
%
% Change log:
%   04 Mar 2013, P.Hehrmann - created
%   06 Aug 2013, PH - added: initial seek of noise; nominal speech RMS now optional
classdef OlsaMixerUnit < ProcUnit
    
    properties
        nomSpeechRms = []; % nominal speech RMS (linear scale)
        lvlType = 'rms'; % string or cell array of strings: 'rms' / 'peak' / 'rel';  
        sensIn = 0;  % input sensitivity (dB SPL equivalent to digital FS)
        wrap = 0; % repeat shorter inputs to match duration of the longest?
        durFade = 0; % duration of cross-fade when wrapping signals [s]
        lvlDb = [0 0]; % 2-element vector of levels per channel in dB
        delays = []; % vector of onset delays for each input [s]
        noiseSeek = 0; % initial seek position within noise [s]; -1 for random
    end
    
    properties(SetAccess=private)
        clip = 0; % flag: did clipping occur during the most recent call of run() ? [0/1]
    end
    
    methods
        function obj = OlsaMixerUnit(parent, ID, nomSpeechRms, lvlType, sensIn, delays, wrap, durFade, noiseSeek)
        %  obj = OlsaMixerUnit(parent, ID, nomSpeechRms, lvlType, sensIn, delays, wrap, durFade)
            
            % declare required shared props
            obj = obj@ProcUnit(parent, ID, 2, 1);
            obj.sharedProps = struct('name', {'sp_fs'}, 'value', {[]});
                     
            assert((isscalar(nomSpeechRms)&&isnumeric(nomSpeechRms)) || isempty(nomSpeechRms),...
                    'nomSpeechRms must be a numeric scalar or empty.');
            obj.nomSpeechRms = nomSpeechRms;
            
            assert(ischar(lvlType), 'lvlType must be a string.');
            obj.lvlType = lvlType;
            
            obj.sensIn = sensIn;
            
            if nargin > 4
                assert(isempty(delays) || (length(delays)==2),...
                    'Length of delays must be 0 or equal to nInputs.');
                assert(all(delays >= 0), 'Delays must be non-negative.');
                obj.delays = delays;
            end
            
            if nargin > 5
                assert(isscalar(wrap), 'wrap must be scalar (logical or 0/1).');
                obj.wrap = wrap;
            end
            
            if nargin > 6
                assert(isscalar(durFade) && durFade >= 0, 'durFade must be a positive scalar.')
                obj.durFade = durFade;
            end
            
            if nargin > 7
                assert(isscalar(durFade) && (durFade >= 0 || durFade == -1) , 'noiseSeek must be a positive scalar or -1.')
                obj.noiseSeek = noiseSeek;
            end
            
            obj.finalizeInit();
        end
        
        function run(obj)
            
            % get all inputs
            X = cell(1,2);
            X{1} = obj.getData('INPUT_1');
            
            if ~obj.getDataUnit('INPUT_2').dataIsEmpty
                X{2} = obj.getData('INPUT_2');
            end
            
            % call processing routine
            [Y clipping] = olsaMixerFunc(X{:}, obj);
            
            % set output
            obj.setData('OUTPUT_1', Y);
            
            % set clipping flag
            obj.clip = clipping;
        end
    end
end