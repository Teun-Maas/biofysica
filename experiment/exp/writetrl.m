function writetrl(fid,Trl)
% WRITETRL(FID,TRL)
%
% Write Trial-line in an exp-file with file identifier FID.
% TRL - Trial Number
%
% See also GENEXPERIMENT, FOPEN, and the documentation of the Auditory
% Toolbox

if nargin<2
	Trl = 0;
end
fprintf(fid,'\n');
fprintf(fid,'%s\n',['% Trial: ' num2str(Trl)]);
fprintf(fid,'%s\n','==>');
