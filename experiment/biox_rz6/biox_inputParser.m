classdef biox_inputParser < inputParser
    % inputParser with slightly modified error handling
    methods
        function this = biox_inputParser
            this = this@inputParser;
        end
        
        function parse(this, varargin)
            try
                parse@inputParser(this, varargin{:});
            catch ME
                switch ME.identifier
                    case 'MATLAB:InputParser:ParamMustBeChar'
                        % The normal error message is 
                        % "Expected a string scalar or character vector for the parameter name."
                        % which is confusing when there are just too many
                        % arguments
                        ME2 = MException( ...
                            'MATLAB:InputParser:ParamMustBeCharOrTooManyArguments', ...
                            'Too many arguments, or parameter name expected.');
                        throw(ME2)
                    otherwise
                        rethrow(ME);
                end
            end
        end
        
    end
    
end