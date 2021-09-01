% ReadWavUnit < ProcUnit
% 
% Read wav data from file. The wav filename can either be specified as
% content of an input DataUnit, or via the property "wavFile" if no input
% unit is defined.
%
% ReadWavUnit Properties:
%   wavFile - wav file name
%
% ReadWavUnit Methods:
%   ReadWavUnit - constructor
%   run - execute processing
%
% Change log:
%   Apr 2012, M.Milczynski - created
%   12 Dec 2012, P.Hehrmann - defaults for nInput, nOutput
%   14 Jan 2012, PH - wav file can now be specified either by an input
%                   DataUnit (previous behavior), or by propterty 'wavFile' 
%                   (new behavior); might deprecate old behavior for
%                   release version
classdef ReadWavUnit < ProcUnit
    
    properties
        wavFile = '';  % wav file name
    end
    methods
        function obj = ReadWavUnit(parent, ID, varargin)
        % Constructor, generate new ReadWavUnit object
        % Option 1: ReadWavUnit(parent, ID [, wavFile])
        % Option 2: ReadWavUnit((parent, ID, nInput, nOutput)) [for backwards compatibility]
            if (nargin < 4) % new syntax: ReadWavUnit(parent, ID [, wavFile ] )
                nInput = 0;
                nOutput = 1;
            else % old syntax for backward compatibility: ReadWavUnit(parent, ID, nInput, nOutput)
                nInput = varargin{1};
                nOutput = varargin{2};                
            end
            
            obj = obj@ProcUnit(parent, ID, nInput, nOutput);
                
            if (nargin == 3)
                obj.wavFile = varargin{1};
            end

            %
            obj.sharedProps(1).name = 'sp_fs';
            obj.sharedProps(1).value = 17400;
            %
            obj.finalizeInit();
        end
        
        
        function run(obj) 
            if obj.inputCount > 0 % input unit is specified (old)
                name = obj.getData('INPUT_1');
                signalIn = readWavFunc(name, obj); 
            else % no input unit is specified (new)
                signalIn = readWavFunc(obj);
            end
            
            obj.setData('OUTPUT_1', signalIn);
        end
    end
end

