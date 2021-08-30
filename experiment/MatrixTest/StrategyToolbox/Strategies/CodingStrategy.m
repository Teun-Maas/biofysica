% CodingStrategy < handle & dynamicprops
% 
% Change log:
%    25 Jun 2012, P.Hehrmann - added GraphChanged event
%    10 Sep 2012, PH - added "verbose" property
%    13 Sep 2012, PH - added "resetNonRootUnits" function
%    29 Jan 2013, M. Milczynski - added name property
%    28 Aug 2013, PH - added call to ProcUnit.updateDepth in connectProcUnits
classdef CodingStrategy < handle & dynamicprops
   properties
       procUnits = {};
       verbose = 0; % 1 to enable text output, 0 to disable
       name = '';
   end
   
   events
       GraphChanged;
   end
   
   methods
       function addProcUnit(obj, unit)
           obj.procUnits{end+1} = unit;
           
           notify(obj,'GraphChanged');
       end
       function connectProcUnits(obj, puIDOut, duIDOut, puIDIn, duIDIn)
           puOut = obj.locateProcUnitByID(puIDOut);
           puIn = obj.locateProcUnitByID(puIDIn);
           puOut.connectToInput(duIDOut, puIn, duIDIn);
           
           puIn.updateDepth(puOut.depth+1);
           notify(obj,'GraphChanged');
       end
       function [pu, i] = locateProcUnitByID(obj, ID)
            found = false;
            for i = 1:length(obj.procUnits)
              if strcmp(obj.procUnits{i}.ID, ID)
                  pu = obj.procUnits{i};
                  found = true;
                  break
              end
            end
            if ~found
               error('no procUnit with ID %s found', ID); 
            end
       end
       function run(obj)
           for i = 1:length(obj.procUnits)
               % get the next procUnit
               puCur = obj.procUnits{i};
               if obj.verbose
                    fprintf(1, 'ProcUnit: %s\n', puCur.ID);
               end
               % execute procUnit
               puCur.run();
               % propagate output to connected units
               puCur.propagateOutput();
           end
       end
	   
	   % Reset the input and output DataUnits of all ProcUnits
       function resetDataUnits(obj)
         for i=1:length(obj.procUnits)
            obj.procUnits{i}.resetDataUnits();
         end
       end       
	   
       % Reset the output DataUnits of all ProcUnits, and the input of all ProcUnits except the roots of the CodingStrategy
	   function resetNonRootUnits(obj)
           % for every ProcUnit in the strategy:
           for k=1:length(obj.procUnits)
               pu = obj.procUnits{k};
               % for every output unit in this ProcUnit:
               for i=1:pu.outputCount
                   du = pu.getDataUnit(sprintf('OUTPUT_%d',i));
                   % clear the output itself
                   du.resetData();
                   % reset all input DataUnits receiving from this output
                   for j=1:du.getNumberConnections()
                       targetInput = du.connection(j, 1).getInputDataUnit();
                       targetInput.resetData();
                   end
               end
           end
           notify(obj,'GraphChanged');
       end % function
       
   end  % methods
   
end % classdef