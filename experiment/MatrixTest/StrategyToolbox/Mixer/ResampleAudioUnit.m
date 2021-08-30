classdef ResampleAudioUnit < ProcUnit
   properties (SetAccess = private)
       fsInVec = [];
       fsOutVec = [];
   end
   methods 
       function obj = ResampleAudioUnit(parent, ID, fsIn, fsOut)
       
           if length(fsIn) ~= length(fsOut)
               error('In- and output vectors with fs specs must have equal length');
           end
           obj = obj@ProcUnit(parent, ID, length(fsIn), length(fsOut));

           % ... specifiy shared properties here
           obj.fsInVec = fsIn;
           obj.fsOutVec = fsOut;
           % [MM]: dirty code!
           if length(fsOut) == 1
               obj.sharedProps = struct('name', {'sp_fs'}, 'value', {fsOut});
           else
               obj.sharedProps = struct('name', {}, 'value', {});
           end
           obj.finalizeInit();
       end
       function run(obj)
           % change this line/add additional commands if necessary                
           for i = 1:obj.inputCount
               slot = num2str(i);
               x = obj.getData(['INPUT_' slot ]);
               if obj.fsInVec(i) ~= obj.fsOutVec(i)
                   x = resample(x, obj.fsOutVec(i), obj.fsInVec(i));
               end
               obj.setData(['OUTPUT_' slot ], x);
           end
       end
   end  
end