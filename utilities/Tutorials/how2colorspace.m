close all
clear all

%% [H,C,L]
H = 260;
C = 0:100;
L = 0:100;
[C,L]	= meshgrid(C,L);
H		= repmat(H,size(C));

[m,n] = size(C);
LCH		= [L(:) C(:) H(:)];
RGB		= pa_LCH2RGB(LCH);
R = reshape(RGB(:,1),m,n);
G = reshape(RGB(:,2),m,n);
B = reshape(RGB(:,3),m,n);
whos R G B RGB
RGB = reshape(RGB,m,n,3);
whos R G B RGB
image(RGB)
axis square;
set(gca,'YDir','normal');
return
%% Sequential colour palettes
ncol	= 64;
l		= linspace(0,1,ncol);
H		= repmat(300,1,ncol);
C		= zeros(1,ncol);
L		= 90-l*30;
LCH		= [L;C;H]';
RGB		= pa_LCH2RGB(LCH);

col		= jet(ncol);
figure(1)
subplot(311)
for ii	= 1:ncol
	plot(ii,1,'ks','MarkerFaceColor',col(ii,:),'MarkerSize',15);
	hold on
	plot(ii,2,'ks','MarkerFaceColor',RGB(ii,:),'MarkerSize',15);
end
% axis square;
xlim([0 ncol+1]);
ylim([0 3])
axis off

%% [H,0+f(i)*Cmax,Lmax-f(i)*(Lmax-Lmin)] 
% with f(i) = i^p is a convenient transformation

ncol	= 64;
l		= linspace(0,1,ncol);
p		= 1;
fl		= l.^p;
Cmax	= 100;
Lmax	= 100;
Lmin	= 30;
H = 150;
H		= repmat(H,1,ncol);
C		= zeros(1,ncol)+fl*Cmax;
L		= Lmax-fl*(Lmax-Lmin);
LCH		= [L;C;H]';
RGB2		= colorspace('RGB<-LCH',LCH);
RGB		= pa_LCH2RGB(LCH);

subplot(312)
for ii	= 1:ncol
	plot(ii,1,'ks','MarkerFaceColor',col(ii,:),'MarkerSize',15);
	hold on
	plot(ii,2,'ks','MarkerFaceColor',RGB(ii,:),'MarkerSize',15);
end
% axis square;
xlim([0 ncol+1]);
ylim([0 3])
axis off


%% [H2-i*(H1-H2),Cmax-f(i)*(Cmax-Cmin),Lmax-f(i)*(Lmax-Lmin)];

ncol	= 64;
l		= linspace(0,1,ncol);
p		= 1;
fl		= l.^p;
Cmax	= 90;
Cmin	= 10; 
Lmax	= 70;
Lmin	= 70;
H1		= 0;
H2		= 360;
H		= H2-l*(H1-H2);
C		= Cmax-fl*(Cmax-Cmin);
L		= Lmax-fl*(Lmax-Lmin);
LCH		= [L;C;H]';
RGB		= colorspace('RGB<-LCH',LCH);
RGB		= pa_LCH2RGB(LCH);
RGB = flipud(RGB);
col		= jet(ncol);

subplot(313)
for ii	= 1:ncol
	plot(ii,1,'ks','MarkerFaceColor',col(ii,:),'MarkerSize',15);
	hold on
	plot(ii,2,'ks','MarkerFaceColor',RGB(ii,:),'MarkerSize',15);
end
% axis square;
xlim([0 ncol+1]);
ylim([0 3])
axis off
