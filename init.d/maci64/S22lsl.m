% Set up matlab path for the labstreaminglayer library
pathdirs = { ...
   [ biofysica_root '/liblsl/maci64' ]
   };

for i = 1:size(pathdirs,1)
   d = pathdirs{i};
   if exist(d,'dir')
      p=genpath(d);
      addpath(p);

      exclude_d=[d '/init.d'];
      exclude_p=genpath(exclude_d);
      rmpath(exclude_p);

      exclude_d=[d '/.svn'];
      exclude_p=genpath(exclude_d);
      rmpath(exclude_p);

      exclude_d=[d '/.git'];
      exclude_p=genpath(exclude_d);
      rmpath(exclude_p);

      exclude_d=[d '/liblsl'];
      exclude_p=genpath(exclude_d);
      rmpath(exclude_p);

      % disp(['added directory ', d, ' to the MATLAB path']);
   else
      disp(['directory ', d, ' does not exist, not added to the MATLAB path']);
   end
   clear pathdirs d p exclude_d i exclude_p
end


