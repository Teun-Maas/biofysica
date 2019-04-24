function retval=biofysica_root
% BIOFYSICA_ROOT is used by the biofysica toolbox initializations scripts
% and returns the root directory of the toolbox.

    persistent loc

    if isempty(loc)
        [loc,~,~]=fileparts(mfilename('fullpath'));
    end
    
    retval=loc;

end

