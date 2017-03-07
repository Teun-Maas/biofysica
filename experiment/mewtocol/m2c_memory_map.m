classdef m2c_memory_map < handle

    properties (Access=protected)

        memory_map;

    end

   
    methods
        function this = m2c_memory_map(filename)
            if verLessThan('matlab', '9.1.0')
                error('m2c_memory_map needs matlab version >= 2016b');
            end
            this.memory_map = readtable(filename);
            % remove rows missing identifier or memory address
            this.memory_map = ...
                rmmissing( this.memory_map, 'DataVariables', { 'Identifier', 'FP_Address' });
            % use 'Identifier' in table as row names
            this.memory_map.Properties.RowNames = this.memory_map.Identifier;
        end
   
        function IEC_struct = lookup(this, Identifier)
            try
                IEC_address = char(this.memory_map{Identifier, 'IEC_Address'});
            catch ME
                IEC_struct = [];
                error('Identifier not found: %s', Identifier);
            end
            IEC_type = char(this.memory_map{Identifier, 'Type'});
            IEC_struct = this.address2struct(Identifier, IEC_address, IEC_type);
        end

        function IEC_struct = address2struct(this, Identifier, IEC_address, IEC_type)
            [fields,n]=sscanf(IEC_address, '%c%c%c%d.%d.%d');
            if (n==5)
              fields=[fields;NaN];
            end
            ts=this.find_matlab_type(IEC_type);

            IEC_struct=struct('location', char(fields(2)), 'type_id', char(fields(3)), ...
                'file', fields(4), 'element', fields(5), 'bit', fields(6), ...
                'offset', 0, 'length', 1, 'identifier', Identifier, ...
                'string', IEC_address, 'type', IEC_type, 'class', ts);
        end

        function ts = find_matlab_type(this,IEC_type)
            switch IEC_type
            case 'BOOL'
                ts='boolean';
            case 'BOOL_OVERLAPPING_DUT'
                ts='uint16';
            case 'INT'
                ts='int16';
            case 'UINT'
                ts='uint16';
            case 'WORD'
                ts='uint16';
            case 'DINT'
                ts='uint32';
            case 'REAL'
                ts='single';
   
            otherwise
                if strncmp(IEC_type, 'ARRAY', 5)
                  ofpos = findstr(IEC_type, 'OF ');
                  tppos = ofpos+3;
                  IEC_type = IEC_type(tppos:end);
                  ts=this.find_matlab_type(IEC_type);
                else
                  ts='unknown';
                end
            end
        end

        function varargout = subsref(this, s)
            switch s(1).type
            case '.'
               if length(s) == 1
                  % Implement this.PropertyName
                  s=this.lookup(s(1).subs);
                  varargout={s};
                  ...
               else
                  varargout = {builtin('subsref',this,s)};
               end
         %  case '()'
         %     if length(s) == 1
         %        % Implement this(indices)
         %        ...
         %     elseif length(s) == 2 && strcmp(s(2).type,'.')
         %        % Implement this(ind).PropertyName
         %        ...
         %     elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
         %        % Implement this(indices).PropertyName(indices)
         %        ...
         %     else
         %        % Use built-in for any other expression
         %        varargout = {builtin('subsref',this,s)};
         %     end
         %  case '{}'
         %     if length(s) == 1
         %        % Implement this{indices}
         %        ...
         %     elseif length(s) == 2 && strcmp(s(2).type,'.')
         %        % Implement this{indices}.PropertyName
         %        ...
         %     else
         %        % Use built-in for any other expression
         %        varargout = {builtin('subsref',this,s)};
         %     end
            otherwise
               error('Not a valid indexing expression')
   end

        end
    end
   
end
