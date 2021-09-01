classdef Strategy < CodingStrategy
    methods
        function obj = Strategy()
            obj = obj@CodingStrategy();
        end
        function update(obj)
            props = properties(obj);
            for i=1:length(props)
                p = obj.findprop([props{i} '_update']);
                if ~isempty(p)
                    propFuncHandle = obj.([props{i} '_update']);
                    obj.(props{i}) = propFuncHandle(obj);
                end
            end
        end
    end
end
