classdef biox_dummy_client < biox_abstract_client

    methods

        function write(this, tagname, value, offset) %#ok<INUSL>
            if nargin < 4
                offset = 0;
            end
            fprintf('biox_dummy_client.write():\n');
            disp(tagname);
            disp(value);
            disp(offset);
            fprintf('end\n');
        end

        function data=read(this, tagname, offset, nWords, nChannels) %#ok<INUSL>
            if nargin < 3
                offset = 0;
            end
            if nargin < 4
                nWords = 1;
            end
            if nargin < 5
                nChannels=1;
            end
            data=[];
            fprintf('biox_dummy_client.read():\n');
            disp(tagname);
            disp(offset);
            disp(nWords);
            disp(nChannels);
            fprintf('end\n');
        end

    end
end
