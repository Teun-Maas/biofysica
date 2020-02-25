classdef rz6_dummy < handle

    methods
        function r=ConnectRZ6(this, str, num)
           disp('ConnectRZ6:');
           disp(varargin); 
           r=true;
        end

        function r=LoadCOF(this,filename)
           disp('LoadCOF:');
           disp(varargin); 
           r=true;
        end

        function r=WriteTagVEX(this, varargin)
           disp('WriteTagVEX:');
           disp(varargin); 
           r=true;
        end

        function r=ReadTagVEX(this, varargin)
           disp('ReadTagVEX:');
           disp(varargin); 
           r=true;
        end

        function r=Run(this)
           disp('Run:');
           r=true;
        end
    end
end
