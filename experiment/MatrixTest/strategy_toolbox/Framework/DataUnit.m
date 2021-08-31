% DataUnit < handle
% 
% Change log:
% 14/09/2012, P.Hehrmann - added function "dataIsEmpty"
classdef DataUnit < handle
    properties (SetAccess = private)
        ID = ''
        type = 'INPUT'; % oder 'OUTPUT'
        connection;
    end
    properties (Access = private)
        connectionCount = 0;
        isConnected = false;
        data = [];
    end
    methods
        function obj = DataUnit(type, ID)
            obj.type = type;
            obj.ID = ID;
        end
        function addConnection(obj, procUnit, id)
            if ~obj.isConnected
                obj.connection = Connection(procUnit, id);
                obj.connectionCount = obj.connectionCount + 1;
                obj.isConnected = true;
            else
                obj.connectionCount = obj.connectionCount + 1;
                obj.connection(obj.connectionCount, 1) = ...
                    Connection(procUnit, id);
            end            
        end
        function data = getData(obj)
           % obj
           if isempty(obj.data)
               error('Data slot empty');
           end
           data = obj.data;
        end
        function setData(obj, data)
           if ~isempty(obj.data)
              error('Trying to set non-empty data slot, reset first!'); 
           end
           obj.data = data;
        end
        function resetData(obj)
            obj.data = [];
        end
        function n = getNumberConnections(obj)
            n = obj.connectionCount; 
        end
        function e = dataIsEmpty(obj)
            e = isempty(obj.data);
        end
    end
    %
    % methods (Abstract)
    %     
    % end
end