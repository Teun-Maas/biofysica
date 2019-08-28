classdef r6z_dummy < handle

    methods
        function ConnectRZ6(this, str, num)
           disp('ConnectRZ6:');
           disp(varargin); 
        end

        function LoadCOF(this,filename)
           disp('LoadCOF:');
           disp(varargin); 
        end

        function WriteTagVEX(this, varargin)
           disp('WriteTagVEX:');
           disp(varargin); 
        end

        function ReadTagVEX(this, varargin)
           disp('ReadTagVEX:');
           disp(varargin); 
        end

        function Run(this)
           disp('Run:');
        end
    end
end
