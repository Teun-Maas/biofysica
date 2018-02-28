% d = dir;
% sel = [d.isdir];
% d = d(sel);
% dnames = {d.name};
% dnames = dnames(3:end);
% 
% 
% ndirs = numel(dnames);
% for ii = 1:ndirs
%     dname = dnames{ii};
%     cd(dname)
%     
%     f = dir('*.m');
%     fnames = {f.name}
%     nfiles = numel(fnames)
%     for jj = 1:nfiles
%            fname = fnames{jj};
%            dest = [fname(1:end-2) 'Minor.m']
% movefile(fname,dest)
%            
%     end
%     cd ..
%     
% end
% 
