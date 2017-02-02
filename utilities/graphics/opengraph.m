function h = openGraph (width,height,mag)
%% Functions for opening and saving graphics that operate the same for
% Windows and Macintosh and Linux operating systems. At least, that's the hope!
%
% See also DBDA2E_UTILITIES, SAVEGRAPH
if nargin<1
	width = 7;
end
if nargin<2
	height = 7;
end
if nargin<3
	mag = 1.0;
end

h = figure;
get(h)


%   if ( .Platform$OS.type != "windows" ) { % Mac OS, Linux
%     tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" ,
%                         ... ) )
%     if ( class(tryInfo)=="try-error" ) {
%       lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
%       graphics.off()
%       X11( width=width*mag , height=height*mag , type="cairo" , ... )
%     }
%   } else { % Windows OS
%     tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
%     if ( class(tryInfo)=="try-error" ) {
%       lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
%       graphics.off()
%       windows( width=width*mag , height=height*mag , ... )
%     }
%   }
% }
%