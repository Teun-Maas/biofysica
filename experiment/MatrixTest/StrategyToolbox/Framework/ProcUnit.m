% ProcUnit < dynamicprops
% 
% Change log:
% 17/09/2012, P.Hehrmann - notify parent strategy of the following graph changes: 
%                          DataUnits being reset, set, or output being propagated
% 06/12/2012, PH - added 'verbose' property (default = 0) to enable/disable
%                  "non-critical" output to the console
% 25/01/2013, PH - bug fix: resetDataUnits
% 28/08/2013, PH - added "depth" property and "updateDepth" method

classdef ProcUnit < dynamicprops
    properties
       parent;
       sharedProps = struct();
       dependentSharedProps = struct('name', {});
       verbose = 0; % provide extended concole output? (1/0)
    end
    properties (SetAccess = private)
        inputCount = 0;
        outputCount = 0;
        ID = '';
        depth = 0;   % depth of this unit within parent Strategy (i.e. longtest possible path to a root)
    end
    properties (Access = private)
        input;
        output;
    end
    properties (Constant = true)
        specialPref = 'sp';
        inputPref = 'INPUT_';
        outputPref = 'OUTPUT_';
    end
    methods
        function obj = ProcUnit(parent, ID, nInputs, nOutputs)
            obj.parent = parent;
            obj.ID = ID;
            if nargin >= 3
                for i=1:nInputs
                    obj.addDataUnit('input'); 
                end
            end
            if nargin == 4
                for i=1:nOutputs
                    obj.addDataUnit('output'); 
                end
            end
            if obj.hasParent()
               parent.addProcUnit(obj); 
            end
        end
    end
    methods (Access = protected)
        function o = hasParent(obj)
            o = ~isempty(obj.parent);
        end
        function updateParent(obj)
            if obj.hasParent()
                obj.parent.update();
            end
        end        
        function finalizeInit(obj)
            obj.registerProps();
            obj.updateParent();
        end 
        function registerShared(obj, propName, propValue, ...
                updateFuncName)
            % Special-property that can only be set by obj.parent
            % - modify name to 'sp_propName' if necessary
            % - add property to obj and obj.parent
            % - make SetAccess private in obj
            % - if this property impacts another property, i.e. in other
            %   property depends on this one, a propUpdateFuncHandle is
            %   provided as an argument and will be executed by the parent
            %   strategy upon change to propValue 
            
            if ~startsWith(propName, obj.specialPref)
               % add sp_ suffix if necessary
               newPropName = [obj.specialPref '_' propName];
            else
               newPropName = propName; 
            end
            % add property to obj and ...
            obj.addprop(newPropName);
            isInStrategy = obj.parent.findprop(newPropName);
            if isempty(isInStrategy)
                % ... to parent if not available yet
                obj.parent.addprop(newPropName);
                obj.parent.(newPropName) = propValue;
                obj.(newPropName) = propValue;
            else
                % otherwise take
                obj.(newPropName) = obj.parent.(newPropName);                
            end
            if nargin == 4 && ~isempty(updateFuncName)
                if iscell(updateFuncName)
                    for i=1:length(updateFuncName)
                        impactedProp = ProcUnit.getImpactedPropName(updateFuncName{i});
                        updateName = [impactedProp '_update'];
                        if isempty(obj.parent.findprop(updateName))
                            obj.parent.addprop(updateName);
                            obj.parent.(updateName) = eval(['@' class(obj) '.' updateFuncName{i}]);
                        end
                    end
                else
                    impactedProp = ProcUnit.getImpactedPropName(updateFuncName);
                    updateName = [impactedProp '_update'];
                    if isempty(obj.parent.findprop(updateName))
                        obj.parent.addprop(updateName);
                        obj.parent.(updateName) = eval(['@' class(obj) '.' updateFuncName]);
                    end
                end
                parentProp = obj.parent.findprop(newPropName);
                parentProp.SetMethod = @setImpacts;
            end
            prop = obj.findprop(newPropName);
            prop.GetMethod = @getAccessSharedSpecial;
            prop.SetAccess = 'private';          
            
            function val = getAccessSharedSpecial(obj)
                % always return parent value
                val = obj.parent.(newPropName);
            end        
            function setImpacts(parent, val)
                if ~isempty(obj.findprop([newPropName '_setMethod']))
                    val = obj.([newPropName '_setMethod'])(val);
                end
                parent.(newPropName) = val;
                parent.update();
            end
        end
        function registerDependentShared(obj, propName)
            % Shared-property that is dependent on other shared properties
            % - modify name to 'sp_propName'
            % - add property to obj and obj.parent
            % - make SetAccess private in both obj and obj.parent
            
            if ~startsWith(propName, obj.specialPref)
               newPropName = [obj.specialPref '_' propName];
            else
               newPropName = propName; 
            end
            obj.addprop(newPropName);
            isInStrategy = obj.parent.findprop(newPropName);
            if isempty(isInStrategy)
                obj.parent.addprop(newPropName);
                obj.parent.(newPropName) = [];
                obj.(newPropName) = [];
                parentProp = obj.parent.findprop(newPropName);
                parentProp.SetAccess = 'private';  
            else
                obj.(newPropName) = obj.parent.(newPropName);                
            end
            prop = obj.findprop(newPropName);
            prop.GetMethod = @getAccessSharedSpecial;
            prop.SetAccess = 'private'; 
            
            function val = getAccessSharedSpecial(obj)
                % always return parent value
                val = obj.parent.(newPropName);
            end        
        end
        function registerProps(obj)
            if ~isempty(fieldnames(obj.sharedProps))
                for i = 1:length(obj.sharedProps)
                    if isfield(obj.sharedProps(i), 'update')
                        obj.registerShared(obj.sharedProps(i).name, ...
                            obj.sharedProps(i).value, ...
                            obj.sharedProps(i).update);
                    else
                        obj.registerShared(obj.sharedProps(i).name, ...
                            obj.sharedProps(i).value);                        
                    end
                end
            end
            if ~isempty(fieldnames(obj.dependentSharedProps))
                for i = 1:length(obj.dependentSharedProps)
                    obj.registerDependentShared(obj.dependentSharedProps(i).name);
                end
            end
        end
    end
    
    methods
        function addDataUnit(obj, type)
            type = lower(type);
            obj.([type 'Count']) = obj.([type 'Count']) + 1;
            id = [obj.([type 'Pref']) num2str(obj.([type 'Count']))];
            if isempty(obj.(type))            
                obj.(type) = DataUnit(upper(type), id);
            else
                obj.(type)(obj.([type 'Count']), 1) = ...
                    DataUnit(upper(type), id);
            end
        end
        
        function connectToInput(obj, objOutID, dest, destInID)
            du = obj.getDataUnit(objOutID);
            du.addConnection(dest, destInID);
        end
        
        function du = getDataUnit(obj, ID)
            tmp = textscan(ID, '%s', 'delimiter', '_');
            type = lower(tmp{1}{1});
            dataUnits = obj.(type);
            found = false;
            for i=1:length(dataUnits)             
                if strcmp(dataUnits(i, 1).ID, ID)
                    du = dataUnits(i, 1);
                    found = true;
                    break;
                end
            end
            if ~found
                error('No dataUnit with ID %s found', ID);
            end
        end
        
        function data = getData(obj, ID)
           du = obj.getDataUnit(ID);
           data = du.getData();
        end
        
        function data = setData(obj, ID, data)
           du = obj.getDataUnit(ID);
           graphChange = du.dataIsEmpty();
           du.setData(data);
           if graphChange 
               notify(obj.parent,'GraphChanged');
           end
        end
        
        function resetDataUnit(obj, ID)
            du = obj.getDataUnit(ID);
            graphChange = ~du.dataIsEmpty();
            du.resetData();
            if graphChange 
               notify(obj.parent,'GraphChanged');
            end
        end
        
        function resetDataUnits(obj)
        % Reset all input and output DataUnits    
            graphChange = false;
            for din = obj.input(:)'
                graphChange = graphChange || ~din.dataIsEmpty();
                din.resetData();
            end
            for dout = obj.output(:)'
                graphChange = graphChange || ~dout.dataIsEmpty();
                dout.resetData();
            end
            if (graphChange)
                notify(obj.parent,'GraphChanged');
            end
        end
        
        function propagateOutput(obj)
            for i=1:obj.outputCount
               du = obj.output(i);
               if strcmp(du.type, 'OUTPUT')
                   data = du.getData();
                   for j=1:du.getNumberConnections()
                       targetInput = du.connection(j, 1).getInputDataUnit();
                       targetInput.setData(data);
                   end
               else
                  error('Cannot propagate data from an input-type data-unit'); 
               end
            end
            notify(obj.parent,'GraphChanged');
        end
        
        function d = updateDepth(obj, minDepth, nPrevSteps)
            % Increase depth of obj to minDepth if necessary and trigger
            % update of all descendants of obj. nPrevSteps is the recursion
            % depth and is used to avoid infinite loops
            if nargin < 3
                nPrevSteps = 0;
            end
            
            if nPrevSteps >= length(obj.parent.procUnits)-1
                error('Maximum recursion depth exceeded: strategy graph must be loopy.');
            end
      
            if minDepth > obj.depth
                obj.depth = minDepth;
                for DU = obj.output(:)'
                    for con = DU.connection(:)'
                       con.procUnit.updateDepth(obj.depth+1, nPrevSteps+1);
                    end
                end
            end
            d = obj.depth;
        end
    end  
    
    methods (Abstract)
        run(obj)
    end
    methods (Static)
        function propName = getImpactedPropName(updateFuncName)
            ptr = strfind(updateFuncName, '_update');
            propName = updateFuncName(1:ptr-1);
        end
    end
end