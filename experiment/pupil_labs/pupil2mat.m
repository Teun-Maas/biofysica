function result=pupil2mat(src,dst)
    % PUPIL2MAT - convert pupil_data file to .mat file using the external
    % 
    %   pupil2mat('/path/to/pupil_data', 'pup.mat') converts the
    %   data in the file /path/to/pupil_data to a .mat-file pup.mat in 
    %   the current working directory.
    %
    % python script pupil2mat.py to be found in the same directory as
    % pupil2mat.m.
    % Needs: python3 installed and in the PATH
    % Needs: python libraries: numpy, scipy and msgpack (use pip3 to install these)
    %
    % example: pupil2mat('~/recordings/2017_11_15/011/pupil_data',...
    %              'pupil_data-2017_11_15-011.mat');
    
    % v1.0 GW/20180517
    [here,~,~]=fileparts(mfilename('fullpath'));
    py_pupil2mat=strcat(here,'/pupil2mat.py');
    if isunix()
        python='python3';
        extrapath=':/usr/local/bin:/opt/bin';
        cmd=['PATH=', getenv('PATH'), extrapath, ' ', python, ' ',...
            py_pupil2mat, ' ', src, ' ', dst];
        result=unix(cmd);
    elseif ispc()
        python='python';
        cmd=[python, ' ', py_pupil2mat, ' ', src, ' ', dst];
        result=system(cmd);
    else
        error('unexpected computer or OS type');
    end
end
