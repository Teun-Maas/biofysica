% AudioMixerUnit < ProcUnit
% Mixes any number of input audio signals into a single output signal. For
% every signal, the target level as well as the target level type have to 
% be specified (abs. rms/abs. peak/rel. to input). sensIn defines the peak 
% SPL equivalent to 0dB FS. wrap controls the behaviour of input signals
% are of unequal length (warp-around or zero-pad to match duration of the 
% primary signal). For wrap = 1, durFade controls the duration of 
% a cosine-shaped cross-fade between the end and beginning.
%
% Properties:
%   sensIn   - input sensitivity (dB SPL equivalent to sinusoid at digital FS)
%   lvlType  - string or cell array of strings: 'rms' / 'peak' / 'rel'; 
%              if string, same type is used for all channels;
%              if cell aray, the size must match n
%   lvlDb    - n-vector of levels per channel in dB; for types 'rms' and
%              'peak', lvlDb(i) is in dB SPL; for type 'rel', lvlDb(i) is
%              in dB relative to the input level.
%   delays   - vector of onset delays for each input [s] 
%   primaryIn - which input determines the length of the mixed output
%               (1..nInputs, or []); if [], the longtest input (including 
%               delay) is chosen as primary.
%   wrap     - repeat shorter inputs to match duration of the longest
%              input? [1/0]
%   durFade  -  duration of cross-fade when wrapping signals [s]
%   clip     - read-only indicator: did clipping occur during the most recent
%              call of run()? [0/1]
%
% See also: audioMixerFunc
%
% Change log:
%   29/08/2012, P.Hehrmann - created
%   09/12/2012, PH - added 'delays' and 'primaryIn' and 'clip'; cf. audioMixerFunc.m
classdef AudioMixerUnit < ProcUnit
    
    properties
        lvlDb; % n-vector of levels per channel in dB
        lvlType = 'rms'; % string or cell array of strings: 'rms' / 'peak' / 'rel';  
        sensIn = 0;  % input sensitivity (dB SPL equivalent to digital FS)
        wrap = 0; % repeat shorter inputs to match duration of the longest?
        durFade = 0; % duration of cross-fade when wrapping signals [s]
        delays = []; % vector of onset delays for each input [s]
        primaryIn = []; % index of primary input ([] = longest)
    end
    
    properties(SetAccess=private)
        clip = 0; % flag: did clipping occur during the most recent call of run() ? [0/1]
    end
    
    methods
        function obj = AudioMixerUnit(parent, ID, nInputs, lvlDb, lvlType, sensIn, delays, primaryIn, wrap, durFade)
            % declare required shared props
            obj = obj@ProcUnit(parent, ID, nInputs, 1);
            obj.sharedProps = struct('name', {'sp_fs'}, ...
                'value', {[]});
            
            % set class properties
            assert(length(lvlDb) == nInputs, 'length(lvlDb) must match nInputs.');
            obj.lvlDb = lvlDb;
            
            if nargin > 4
                assert(ischar(lvlType) || ...
                      (iscell(lvlType) && length(lvlType) == nInputs) && all(cellfun(@ischar,lvlType)), ...
                       'lvlType must be either a string of a cell array of strings of length nInputs.');
                obj.lvlType = lvlType;
            end
            
            if nargin > 5
                assert(isscalar(sensIn), 'sensIn must be a scalar.');
                obj.sensIn = sensIn;
            end
            
            if nargin > 6
                assert(isempty(delays) || (length(delays)==nInputs),...
                    'Length of delays must be 0 or equal to nInputs.');
                assert(all(delays >= 0), 'Delays must be non-negative.');
                obj.delays = delays;
            end
            
            if nargin > 7
                assert(isempty(primaryIn) || (primaryIn > 0 && primaryIn <= nInputs && ~mod(primaryIn,1)), ...
                    'primaryIn must be empty or an integer between 1 and nInputs.');
                obj.primaryIn = primaryIn;
            end        
            
            if nargin > 8
                assert(isscalar(wrap), 'wrap must be scalar (should be logical or 0/1).');
                obj.wrap = wrap;
            end
            
            if nargin > 9
                assert(isscalar(durFade) && durFade >= 0, 'durFade must be a positive scalar.')
                obj.durFade = durFade;
            end
            
            obj.finalizeInit();
        end
        
        function run(obj)
            
            % get all inputs
            X = cell(1,obj.inputCount);
            for iIn = 1:obj.inputCount
                tmp = obj.getData(sprintf('INPUT_%d',iIn));
                tmp = tmp(:);
                X{iIn} = tmp; 
            end
            
            % call processing routine
            [Y clipping] = audioMixerFunc(X{:}, obj);
            
            % set output
            obj.setData('OUTPUT_1', Y);
            
            % set clipping flag
            obj.clip = clipping;
        end
    end
end