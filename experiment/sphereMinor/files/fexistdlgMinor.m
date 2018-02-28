function cfg = fexistdlg(cfg)
% [FNAME,DNAME] = FEXISTDLG(FNAME,DNAME)
%
% Check whether the current data file name does not already exist. If so,
% prompts experimenter to continue, stop, or rename.
%
% See also GETFNAME

% 2015 Marc van Wanrooij
% e: marcvanwanrooij@gmail.com
fileexist		= exist([cfg.dname filesep cfg.fname],'file');
while fileexist
	choice1 = 'Continue/overwrite';
	choice2 = 'Stop';
	choice		= questdlg(['File ' cfg.fname ' already exists. What dou you want to do?'],'Check filename',choice1,choice2,choice2);
	switch choice
		case(choice1)
			disp('Continuing...');
			fileexist = false;
		case(choice2)
			% 			disp('Quitting...');
			cfg.fname = [];
			cfg.dname = [];
			return
	end
end
