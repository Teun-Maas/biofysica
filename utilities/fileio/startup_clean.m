%
% STARTUP
%

%disp('running startup.m');

pathdirs = { ...
   'd:/toolbox/PsychToolbox'
   '/Users/gunter/Projects/DPX'
   '/Users/gunter/Projects/Tactile Hand'
   '/Users/gunter/Documents/Matlab/PupilLabs/cosy-zeromq-matlab/zmq'
   '/Users/gunter/Documents/Matlab/PupilLabs/msgpack-matlab'
   '/Users/gunter/Documents/Matlab/PupilLabs/jsonlab'
   };

for i = 1:size(pathdirs,1)
   d = pathdirs{i};
   if exist(d,'dir')
      p=genpath(d);
      addpath(p);

      exclude_d=[d '/.svn'];
      exclude_p=genpath(exclude_d);
      rmpath(exclude_p);

      disp(['startup: added directory ', d, ' to the MATLAB path']);
   else
      disp(['startup: directory ', d, ' does not exist, not added to the MATLAB path']);
   end
end


