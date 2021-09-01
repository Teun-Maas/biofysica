classdef SplitNResampleUnit < ProcUnit
    properties (SetAccess = private)
        splitResVec = [];
    end
    methods
        function obj = SplitNResampleUnit(parent, ID, splitResVec)
            % splitResMatrix has N entries each with a particular sampling
            % frequency
            obj = obj@ProcUnit(parent, ID, 1, length(splitResVec));
            obj.splitResVec = splitResVec;
            obj.sharedProps = struct('name', {'fs'}, 'value', {17400});
            obj.finalizeInit();
        end
        function run(obj)
            
            x = obj.getData('INPUT_1');
            for i=1:obj.outputCount
                if obj.splitResVec(i) ~= obj.sp_fs
                    r = resample(x, obj.splitResVec(i), obj.sp_fs);
                else
                    r = x; 
                end
                obj.setData(['OUTPUT_' num2str(i)], r);
            end
        end
    end
end