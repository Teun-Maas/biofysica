classdef ResampleUnit < ProcUnit
   properties (SetAccess = private)
      % note there is no pu_ prefix 
      resRatio 
      startIdx
   end
   methods
       function obj = ResampleUnit(parent, ID, nInputs, nOutputs, ...
               resRatio, startIdx)
           obj = obj@ProcUnit(parent, ID, nInputs, nOutputs);
           obj.resRatio = resRatio;
           if nargin == 5
              obj.startIdx = resRatio;
           elseif nargin == 6
               obj.startIdx = startIdx;
           end
       end
       function run(obj)
          x = obj.getData('INPUT_1');
          y = resampleFunc(x, obj);
          obj.setData('OUTPUT_1', y);
       end
   end
end