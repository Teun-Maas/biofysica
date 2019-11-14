function update_biofysica
%UPDATE_BIOFYSICA -- update the biofysica toolbox from the gitlab repository  
%You need a working git command in your system PATH for this

   oldcwd = cd;
   root = biofysica_root();
   cd(root);
  
   if strcmp(cd,root)
       fprintf('I found your biofysica toolbox is in ''%s''\n',root);
       input('press <enter> to continue, or <control-C> to interrupt:');
       fprintf('executing git pull to update your biofysica toolbox from the gitlab repository\n');
       status=system('git pull');
       cd(oldcwd);
       if status ~= 0
           error('oops, apparently something went wrong....');
       end
   else
       error('could not change directory to biofysica root ''%s''',root);
   end
   fprintf('done\n');
end