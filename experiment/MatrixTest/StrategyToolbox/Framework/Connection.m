classdef Connection < handle
    properties (SetAccess = private)
       procUnit;
       inputID;
    end
    methods 
        function obj = Connection(pu, id)
            obj.procUnit = pu;
            obj.inputID = id;
        end
        function du = getInputDataUnit(obj)
            du = obj.procUnit.getDataUnit(obj.inputID);
        end
    end
end