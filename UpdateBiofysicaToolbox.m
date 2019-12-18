function UpdateBiofysicaToolbox
%UPDATEBIOFYSICATOOLBOX -- update the biofysica toolbox from the gitlab repository  
%You need a working git command in your system PATH for this

   oldcwd = cd;
   try
       root = biofysica_root();
   catch
       % actually, we should never get here...
       error('Sorry, I can''t find your biofysica toolbox...');
   end
   cd(root);
  
   if strcmp(cd,root)
       fprintf('I found your biofysica toolbox is ''%s''\n',root);
       input('Press <enter> to update, or <control-C> to interrupt:');
       fprintf('\nExecuting ''git pull'' to update your biofysica toolbox from the gitlab repository\n');
       status=system('git pull');
       cd(oldcwd);
       if status ~= 0
           error('oops, apparently something went wrong....');
       end
   else
       error('could not change directory to biofysica root ''%s''',root);
   end
   fprintf('Done.\nYou may need to restart Matlab.\n');
end
