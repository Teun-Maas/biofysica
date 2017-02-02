function twinfit_example
close all

Subject		=	'Tyr/';
% for m = 1:Ncue
% 	for k = 1:Nfmod
% 		switch Ucue(m)
[RT,SOA] = getdata(Subject);


% keyboard
rt	= [RT(1,1,:).mu];


soa = squeeze(SOA(1,1,:));
% some data manipulation
% rt([4 6]) = rt([4 6])-12;
% rt([3 7]) = rt([3 7])+10;
% rt([2 8]) = rt([2 8])-12;
rt1 = RT(1,1,end).rt;
rt2 = RT(1,2,end).rt; % visual

%%
close all
subplot(121)
plot(soa,rt,'ko','MarkerFaceColor','w')

subplot(122)
[f,xi] = ksdensity(rt1);
plot(xi,f);

hold on

[f,xi] = ksdensity(rt2);
plot(xi,f);


%%


% axis square;
% box off
% set(gca,'TickDir','out');
% verline;
%% selection

sel = rt1>150;
rt1 = rt1(sel);

sel = rt2>150;
rt2 = rt2(sel);

%%
rt = rt/1000;
rt1 = rt1/1000;
rt2 = rt2/1000;
samples = twinfit(soa,rt,rt1,rt2,'showPost',true,'showPred',true,'showDiag',true);

function [RT,SOA] = getdata(Subject)
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

%% Data are stored in a matrix:
%
%  1   2     3    4     5	  6		7	  8	   9  10   11	   12		13	   14 15	16
% fmod dens block trl tswitch visint type delay cue RT correct RTstart moddepth day set subject
D		=	getdat(Name);

[RT,SOA] = getRT(D,VisInt,ModDep,Fmod);


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

function [P,Crit]	= getpercentile(D,Per)

if( nargin < 2 )
	Per	=	[25 75];
end

if( isempty(D) )
	P	=	nan(2,1);
	Crit=	NaN;
	return
end

sd		=	sort( D );
P(1,1)	=	round( (Per(1)/100) * length(D) + .5 );
P(2,1)	=	round( (Per(2)/100) * length(D) + .5 );
P(1,1)	=	sd(P(1),1);
P(2,1)	=	sd(P(2),1);

IQR		=	P(2) - P(1);
Crit	=	IQR * 1.5;	%-- Tuckey criterion for outliers --%

function [RT,SOA] = getRT(Dat,VisInt,MD,Fmod)



%  1   2     3    4     5	  6		7	  8	   9  10   11	   12		13	   14 15	16		--%
% fmod dens block trl tswitch visint type delay cue RT correct RTstart moddepth day set subject	--%

%% All A, V and AV data
% colum 7, 0 = auditory, 1 = visual, 2 = audiovisual
selrt = Dat(:,12)>0;
selav		=	Dat(:,7) == 2 & Dat(:,6) == VisInt;


AuVi	=	Dat(selav&selrt,:);


%%
Udelay	=	unique(Dat(:,8));
Ndelay	=	length(Udelay);

% The task cue, 0 = , 1 = , 2 =
Ucue	=	unique(Dat(selav,9));
Ncue	=	length(Ucue);

% Frequency modulation
Nfmod	=	length(Fmod);

%% graphics



cnt		=	1;

%% Unimodal visual
%  1   2     3    4     5	  6		7	  8	   9  10   11	   12		13	   14 15	16
% fmod dens block trl tswitch visint type delay cue RT correct RTstart moddepth day set subject

RT(Ncue,Nfmod,Ndelay).mu	= NaN;
SOA	= NaN(Ncue,Nfmod,Ndelay);
for m = 1:Ncue
	for k = 1:Nfmod
		switch Ucue(m)
			case 0 % auditory target
				% 			ec		=	'r';
				SOA(m,k,:)	=	Udelay;
			case 1 % visual target
				SOA(m,k,:)	=	-1*Udelay;
			case 2 % redundant target
				SOA(m,k,:)	=	Udelay ;
		end
		
		for n = 1:Ndelay
			sel			=	AuVi(:,13) == MD & AuVi(:,1) == Fmod(k) & AuVi(:,9) == Ucue(m) & AuVi(:,8) == Udelay(n) & ~isnan(AuVi(:,10));
			RT(m,k,n).mu	= nanmean(AuVi(sel,12));
			RT(m,k,n).rt	= AuVi(sel,12);
		end
	
	end
	
	cnt = cnt+1;
end

