function pa_spk_checktrial(fname,datname)
% PA_SPK_CHECKTRIAL(FNAME,DATNAME)
%
% CHECK TRIAL ORDER
%
% BrainWare saves data according to stimulus parameters, not to actual
% trial order. The reaction time data however, is stored by actual trial
% order. This function checks whether the trial order as determined in
% PA_SPK_READSRC is correct.
%
% See also PA_SPK_READSRC, PA_SPK_READDAT
%

% 2013 Marc van Wanrooij
% e-mail:marcvanwanrooij@neural-code.com

if nargin<1
	cd('E:\DATA\Test\check order');
	fname	= 'joe6702c01b02.src';
	datname = 'joe6702.dat';
end
[spk,cfg]						= pa_spk_readsrc(fname); %#ok<*NASGU>
[RT, TN, D, A, RV, RD, MD, RW]	= pa_spk_readdat(datname); %#ok<NASGU,ASGLU>

%%
ord = [spk.trialorder];

% check whether ripple parameters are the same
sv = [spk(ord).stimvalues];
rv = sv(5,:)';

rd = round(sv(6,:)'*10)/10;
RD = round(RD*10)/10;
disp([RV rv RD rd]);

 a = [any(RV~=rv) any(RD~=rd)];
 disp(a);

