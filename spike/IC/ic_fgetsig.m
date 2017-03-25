function [ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10]=ic_fgetsig(DATfiles,Nchan,Nsample,ChNr,TrialNr,On,Off);

% [ch1,...ch10]=fgetsig(DATfiles,Nchan,Nsample,ChNr,TrialNr,On,Off)
%   Read signal blocks (Nsample in length) from HUM,AAP,HV or HVN 
%   files in trials TrialNr from sample number On to Off. If On 
%   and/or Off are scalars the on and/or offset of the blocks are 
%   equal in each trial. If the block sizes are unequal, NaN's 
%   fill up the matrix.
%
%   Jeroen Goossens

% number of files to read from
Nfile = size(DATfiles,1);
% number of blocks per file
Nblock = length(TrialNr);

% set on and offsets of blocks
if length(On)==1,
   On  = ones(Nblock,1)*On;
end; 
if length(Off)==1,
   Off = ones(Nblock,1)*Off;
end;       
Npnt = max(Off-On+1);
Indx = 1;

%%%%% init ch matrices
for i=1:length(ChNr),
    cmd = [ 'ch' num2str(i) ' = NaN*ones(Nblock*Nfile,Npnt);' ];
    eval(cmd);
end;

for f=1:size(DATfiles,1),

   %%%%% open data file
   DATfile = DATfiles(f,:);
   fid = fopen(DATfile,'r','l');
   if fid==-1,
     disp(['   Error reading ' DATfile]);
     return;
   end;
   frewind(fid);

   % find file type
   Ext = DATfile(findstr(DATfile,'.')+1:length(DATfile));
   Ext = upper(Ext);
   Nbyte=0;
   if strcmp(Ext,'HUM'), Nbyte=2; end;
   if strcmp(Ext,'AAP'), Nbyte=2; end;
   if strcmp(Ext,'HV'),  Nbyte=4; end;
   if strcmp(Ext,'HVN'), Nbyte=4; end;
   if Nbyte==0, 
     disp(['   Error : Invalid DATfile extension ',Ext]);
     return;
   end;   


   %%%%% read blocks

   for i=1:Nblock,
     % set file pointer at beginning of block
     fpos   = Nbyte*Nchan*(Nsample*(TrialNr(i)-1)+(On(i)-1));
     status = fseek(fid,fpos,'bof');
     if status ~= 0 ,
       disp(['   Error : Invalid file position ' DATfile ' : Trial #' num2str(TrialNr(i))]);
       return;
     end;

     % read block
     Npnt    = Off(i)-On(i)+1;
     if Nbyte==2,
       [mtx,n] = fread(fid,[Nchan,Npnt],'ushort');
     end;
     if Nbyte==4,
       [mtx,n] = fread(fid,[Nchan,Npnt],'float');
     end;
     if n~=Npnt*Nchan,
       disp([' Waring : Trial #' num2str(TrialNr(i)) ' in file ' DATfile ' truncated']);
     end;

     % extract channels
     ns = size(mtx,2);
     for j=1:length(ChNr),
       nr  = ChNr(j);
       cmd = [ 'ch' num2str(j) '(Indx,1:ns) = mtx(nr,:); ' ];
       eval(cmd);
     end;
     Indx = Indx +1;

   end;

   % close file
   fclose(fid);
end;
