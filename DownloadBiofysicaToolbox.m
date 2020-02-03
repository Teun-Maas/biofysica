function DownloadBiofysicaToolbox(targetdirectory)
% DOWNLOADBIOFYSICATOOLBOX - Download and install a copy of the Biofysica Toolbox
% 
%  DownloadBiofysicaToolbox - install in the default location, on MacOs and Windows
%  systems this is in your personal MATLAB directory Documents->MATLAB->biofysica
%  On other UNIX systems this may be in ~/MATLAB/biofysica
%
%  SEE ALSO: UpdateBiofysicaToolbox
    
% GW/20191218 - initial version

    repository='https://gitlab.science.ru.nl/marcw/biofysica.git';
    
    if nargin < 1
        targetdirectory=fullfile(get_matlabdir(),'biofysica');
    else
        targetdirectory=fullfile(targetdirectory);
    end
    
    
    % Find out if a copy of the toolbox is already installed
    btinstalled = true;
    try
        root = biofysica_root();
    catch
        btinstalled = false;
    end
    
    if btinstalled
        fprintf('\nHmmm. You already have a copy of the biofysica toolbox installed\n');
        fprintf('in this directory: %s\n', root);
        fprintf('That old toolbox should be removed before we install a new one\n\n');
        fprintf('To remove the old toolbox delete the directory in %s, and remove the\n', root);
        fprintf('the line with ''addpath'' + ''init_biofysica'' in %s\n\n', which('startup.m'));
        fprintf('If you want to keep the old installation, just remove the line with\n');
        fprintf('''addpath'' and ''init_biofysica'' from %s\n\n',which('startup.m'));
        fprintf('Instead of installing you may also consider to update using the ''UpdateBiofysicaToolbox'' command\n\n');
        fprintf('Then RESTART matlab.\n\n');
        %TEST        error('Old toolbox is installed. Please remove or update it.');
    end
    
    if exist(targetdirectory,'dir') && ~isemptydir(targetdirectory)
        error( ...
            'destination path ''%s'' already exists and is not an empty directory', ...
            targetdirectory);
    end
    
    
    % prerequisites
    % 1. Git
    if ispc
        [status,result] = system('where git');
    else % UNIX
        [status,result] = system('which git');
    end
    if status
        fprintf('To download the biofysica toolbox you need a working git client\n');
        fprintf('The ''git'' command is not in not in your system PATH.\n');
        fprintf('Please download and install a recent git client.\n');
        if ispc
            fprintf('For Windows download and install the git client from\n');
            fprintf('https://git-scm.com/download/win and then run %s again.\n', mfilename);
        elseif ismac
            fprintf('For MacOS download and install the git client from\n');
            fprintf('https://git-scm.com/download/mac and then run %s again.\n', mfilename);
        else % guess this is Unix
            fprintf('For Unix/Linux install the git client using your package manager\n');
        end
        error('Git client is missing. Please install it.');
        
    end
   
    gitcmd = strip(result);
    if ispc
        gitcmd = ['"' gitcmd '"'];
    end
    
    % prerequisites
    % 2. Python 3
    [v, ~, ~] = pyversion;
    fv = str2double(v);
    if fv < 3.8
        fprintf('We need to use a python version >= 3.7\n');
        fprintf('Matlab currently uses:\n');
        pyversion;
        fprintf('\nFinding python versions installed on your system...\n\n');
        if ispc
            [r,str]=system('where /R "C:\Program Files" python');
        else
            [r,str]=system('bash -c ''type -ap python3 python''');
        end
        if r
            fprintf('No python found on your system\n');
        else
            if isspace(str(end))
                str=str(1:end-1);
            end
            str=strsplit(str,'\n');
            for p=str
                [~,ver]=system(sprintf('%s -V',p{1}));
                if isspace(ver(end)), ver=ver(1:end-1); end
                
                fprintf('%s is %s\n',p{1}, ver);
                
                % ver is a string like 'Python 3.7.1'
                ver=strsplit(ver); % get rid of the 'Python' part
                ver=ver{2};
                iver=sscanf(ver,'%i.%i'); % extract major and minor version number
                if iver(1) == 3 && iver(2) >= 7
                    fprintf('To change the Python version to %s, restart Matlab and type this command:\n', ver);
                    fprintf('  pyversion %i.%i\n\n', iver(1), iver(2));
                    fprintf('Then run %s again.\n', mfilename);
                    fprintf('(You will only have to do this once)\n');
                    %TEST                    return
                    break
                end
            end
        end
        fprintf('Please install a python version >= 3.7\n\n');
        fprintf('I suggest to do this first, but you may continue installing without Python now.\n');
        fprintf('However,you may run into trouble later.\n');
        reply = input('Type ''continue'' if you want to continue anyway: ','s');
        if ~strcmp(reply,'continue')
            return;
        end
        fprintf('Proceeding anyway...\n');
    end
    
    fprintf('\nNow installing biofysica toolbox into %s\n', targetdirectory);
    fprintf('Use your science.ru.nl login credentials to start the download\n');
    
    try
        cmd=[gitcmd, ' clone ', repository, ' ', targetdirectory];
        system(cmd);
        %TEST mkdir(targetdirectory);
    catch
        error('Oops, something went wrong here....');
    end
    
    p=what(targetdirectory);
    fprintf('\n\nThe biofysica toolbox has been installed.\n');
    fprintf('Now please check the file ''startup.m'' in %s\nand add or change the following two lines:\n\n',get_matlabdir());
    fprintf('  addpath(''%s'');\n',p.path);
    fprintf('  init_biofysica;\n');
    fprintf('\n');
    fprintf('Done.\n');  
end

function m = get_matlabdir
    if ispc || ismac
        m = fullfile(get_homedir,'Documents','MATLAB');
    else
        m = fullfile(get_homedir,'MATLAB');
    end
    if ~exist(m,'dir')
        warning('Oops, cannot find Matlab userdir...');
        m=[];
    end
end

function h = get_homedir
    if ispc
        h=strcat(getenv('HOMEDRIVE'),'\',getenv('HOMEPATH'));
    else
        h=getenv('HOME');
    end
end


function e = isemptydir(path)
    e = isempty(dir(path));
end


