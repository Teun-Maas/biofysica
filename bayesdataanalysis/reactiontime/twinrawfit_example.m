function twinrawfit_example
close all

Subject		=	'Peter/';
% for m = 1:Ncue
% 	for k = 1:Nfmod
% 		switch Ucue(m)
[RT,MF,Cue,SOA] = getdata(Subject);

%%
% 	Fmod		=	2.^([2; 4; 6]);
	Fmod		=	2.^([1; 4; 8]);

uMF = Fmod;
nMF = numel(uMF);
uSOA = unique(SOA);
nSOA = numel(uSOA);
selcue = Cue==0;
cnt = 0;

nMF
nSOA
for ii = 1:nMF
	for jj = 1:nSOA
		cnt = cnt+1;
		sel = selcue & SOA==uSOA(jj) & MF==uMF(ii);
		
		if sum(sel)
subplot(1,3,ii)
[h,f,xi] = bf_hist(RT(sel),'xi',0:700,'function','pdf','yoffset',jj*0.1);
% [M(ii,jj),I(ii,jj)] = max(f);
% X(ii,jj) = xi(I(ii,jj));
hold on
str = ['MF = ' num2str(uMF(ii))];
title(str)

str = ['MF = ' num2str(uMF(ii))];
text(0,jj*0.002,num2str(uSOA(jj)));
		end
	end
end


% keyboard
return
%%
samples = twinfit(soa,rt,rt1,rt2,'showPost',true,'showPred',true,'showDiag',true);

function [RT,MF,Cue,SOA] = getdata(Subject)
%-- Preprocessing of handlebar data																	--%
%-- Data are stored in a matrix:																	--%
%--
%--	  1   2     3    4     5	  6		7	  8	   9  10   11	   12		13	   14 15	16		--%
%-- fmod dens block trl tswitch visint type delay cue RT correct RTstart moddepth day set subject	--%

% clear all
close all
% clc
tic

%-- Flags --%



%-- Variables --%
Pname		=	'/Users/marcw/DATA/AVAM/';
% Subject		=	'Human3/';
% Subject		=	'Yoolla/';
% Subject		=	'Peter/';

VisInt		=	2.5;
ModDep		=	.25;

if( strcmpi(Subject,'Tyr/') )
	Fmod		=	2.^([2; 4; 6]);
	ExpDate		=	[{'20_07_2015/'};{'21_07_2015/'};{'23_07_2015/'}; ...
		{'27_07_2015/'};{'28_07_2015/'};{'29_07_2015/'}; ...
		{'30_07_2015/'};{'31_07_2015/'};{'03_08_2015/'}; ...
		{'04_08_2015/'};{'05_08_2015/'};{'06_08_2015/'}; ...
		{'07_08_2015/'};{'10_08_2015/'};{'11_08_2015/'}; ...
		{'12_08_2015/'};{'13_08_2015/'};{'14_08_2015/'}; ...
		{'17_08_2015/'};{'18_08_2015/'};{'19_08_2015/'}; ...
		{'20_08_2015/'};{'21_08_2015/'};{'24_08_2015/'}; ...
		{'25_08_2015/'};{'26_08_2015/'};{'27_08_2015/'}];		%-- cell array of strings or []; --%
elseif( strcmpi(Subject,'Peter/') )
	Fmod		=	2.^([1; 4; 8]);
	ExpDate		=	[{'12_08_2015/'};{'13_08_2015/'};{'14_08_2015/'}; ...
		{'17_08_2015/'};{'18_08_2015/'};{'19_08_2015/'}; ...
		{'20_08_2015/'};{'21_08_2015/'};{'24_08_2015/'}; ...
		{'25_08_2015/'};{'26_08_2015/'};{'27_08_2015/'}];
elseif( strcmpi(Subject,'Yoolla/') )
	Fmod		=	2.^([1; 4; 8]);
	ExpDate		=	[{'13_08_2015/'};{'14_08_2015/'};{'17_08_2015/'}; ...
		{'18_08_2015/'};{'19_08_2015/'};{'20_08_2015/'}; ...
		{'21_08_2015/'};{'25_08_2015/'};{'26_08_2015/'}; ...
		{'27_08_2015/'}];
elseif( strcmpi(Subject,'Human3/') )
	Fmod		=	2.^([1; 4; 8]);
	ExpDate		=	[];
end

Aname	=	[Pname Subject];

if( isempty(ExpDate) )
	Name	=	getfolders(Aname);
else
	Nname	=	length(ExpDate);
	Name	=	cell(Nname,1);
	for k=1:length(ExpDate)
		Name(k,1)	=	{[Aname ExpDate{k}]};
	end
end


D			=	getdat(Name);
[RT,MF,Cue,SOA] = getRT(D,VisInt);


function Pname		= getfolders(Name)

Fname	=	dir([Name '*_*']);
Nname	=	length(Fname);

sel		=	nan(Nname,1);
Pname	=	cell(Nname,1);

for k=1:Nname
	if( ~Fname(k,1).isdir )
		sel(k,1)	=	1;
	else
		sel(k,1)	=	0;
	end
	
	if( isunix )
		Pname(k,1)		=	{[Name Fname(k,1).name '/']};
	else
		Pname(k,1)		=	{[Name Fname(k,1).name '\']};
	end
end

sel			=	logical(sel);
Pname(sel)	=	[];

DateNum		=	getdate(Pname);

[~,idx]		=	sortrows(DateNum);
Pname		=	Pname(idx,1);

function DateNum	= getdate(Pname)

Ndat	=	size(Pname,1);
DateNum	=	nan(Ndat,3);
for k=1:Ndat
	P	=	Pname{k,1};
	
	if( isunix )
		idx	=	strfind(P,'/');
	else
		idx	=	strfind(P,'\');
	end
	
	DMY	=	P(idx(end-1)+1:idx(end)-1);
	
	idx	=	strfind(DMY,'_');
	D	=	str2double( DMY(1:idx(1)-1) );
	M	=	str2double( DMY(idx(1)+1:idx(2)-1) );
	Y	=	str2double( DMY(idx(2)+1:end) );
	
	DateNum(k,:)	=	[Y,M,D];
end

function D			= getdat(Name)

Nname	=	size(Name,1);
D		=	[];
for k=1:Nname
	N	=	[Name{k,1} 'Processed' filesep];
	
	F	=	dir([N '*.mat']);
	load([N F.name])
	
	if( size(Mtx,2) == 15 )
		Mtx	=	[Mtx(:,1:8) Mtx(:,7) Mtx(:,9:15)];
	end
	
	D	=	[D; Mtx];											%#ok<AGROW>
end


% Note: datasets 87 & 88 are not good. There was a problem with the LED --%
% Note: dataset 92 is not good. First orange LED dataset; not motivated --%
sel		=	D(:,14) > 0;
D		=	D(sel,:);

% Date 26/05/2015 or set 86: made a type-o '600600' instead of '600 600'--%
sel			=	D(:,8) > 1000;
D(sel,8)	=	ones(sum(sel),1) * 600;
sel			=	D(:,8) < -1000;
D(sel,8)	=	ones(sum(sel),1) * -600;

sel			=	D(:,8) == 74;
D(sel,8)	=	ones(sum(sel),1) * 75;
sel			=	D(:,8) == -74;
D(sel,8)	=	ones(sum(sel),1) * -75;


function [RT,MF,Cue,SOA] = getRT(Dat,VisInt)



%  1   2     3    4     5	  6		7	  8	   9  10   11	   12		13	   14 15	16		--%
% fmod dens block trl tswitch visint type delay cue RT correct RTstart moddepth day set subject	--%

%% All A, V and AV data
% colum 7, 0 = auditory, 1 = visual, 2 = audiovisual
selrt	= Dat(:,12)>100;
selav	=	Dat(:,7) == 2 & Dat(:,6) == VisInt & ~isnan(Dat(:,10));
AV		=	Dat(selav&selrt,:);
RT		= AV(:,12);
MF		= AV(:,1);
Cue		= AV(:,9);
SOA		= AV(:,8);
