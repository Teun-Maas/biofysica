function result=pldata2mat(datadir, topic ,matfile)
    % PLDATA2MAT - convert pupil_data file to .mat file using the external
    % 
    %   pldata2mat('/path/to/pupil_data', TOPIC, MATFILE) converts the
    %   data on TOPIC in the directory /path/to/pupil_data to a .mat-file pup.mat in 
    %   the current working directory.
    %   TOPIC can be 'gaze', 'notify' or 'pupil'.
    %
    % python script pldata2mat.py to be found in the same directory as
    % pldata2mat.m.
    % Needs: python3 installed and in the PATH
    % Needs: python libraries: numpy, scipy and msgpack (use pip3 to install these)
    %
    % example: pldata2mat('~/recordings/2017_11_15/011/',...
    %              'pupil', 'pupil_data-2017_11_15-011.mat');
    
    % v1.0 GW/20190704
    [here,~,~]=fileparts(mfilename('fullpath'));
    py_pldata2mat=strcat(here,'/pldata2mat.py');
    if isunix()
        python='python3';
        %extrapath=':/usr/local/bin:/opt/bin';
        %cmd=['PATH=', getenv('PATH'), extrapath, ' ', python, ' ',...
        %    py_pldata2mat, ' ', datadir, ' ', topic, ' ', matfile];
        extrapath='/usr/local/bin:/opt/bin:';
        cmd=['PATH=', extrapath, getenv('PATH'), ' ', python, ' ',...
            py_pldata2mat, ' ', datadir, ' ', topic, ' ', matfile];
        result=unix(cmd);
    elseif ispc()
        python='python';
        cmd=[python, ' ', py_pldata2mat, ' ' datadir, ' ', topic, ' ', matfile];
        result=system(cmd);
    else
        error('unexpected computer or OS type');
    end
end
