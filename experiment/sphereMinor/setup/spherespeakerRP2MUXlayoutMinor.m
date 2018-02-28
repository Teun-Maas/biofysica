close all;
cfg = spherelookup;

cfg.lookuplabel

RP2 = cfg.lookup(:,2);
uRP2 = unique(RP2);
nRP2 = numel(uRP2);
MX			= cfg.lookup(:,3);
sel			= ismember(MX,[1 2]);
MX(sel)		= 1;
MX(~sel)	= 2;
uMUX		= unique(MX);
nMUX		= numel(uMUX);

AZ		= cfg.lookup(:,5);
EL		= cfg.lookup(:,6);

		
col = jet(nMUX*nRP2);
figure(1)
subplot(121)
hold on
cnt = 0;
for ii = 1:nRP2
	for jj = 1:nMUX
		cnt = cnt+1;
		
		sel		= MX == uMUX(jj) & RP2==uRP2(ii);
		az		= AZ(sel);
		el		= EL(sel);
		chn		= cfg.lookup(sel,1);
		rp2		= RP2(sel);
		mx		= MX(sel);
		
		rpstr = num2str(unique(rp2));
		mxstr = num2str(unique(mx));
		str = ['R' rpstr ', M' mxstr];
		plot(az,el,'ko','MarkerFaceColor',col(cnt,:),'MarkerSize',10);
% 		text(az,el,str,'HorizontalAlignment','center');
	end
end
xlim([-120 120]);
ylim([-120 120]);
axis square
box off
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');


figure(1)
subplot(122)
hold on
cnt = 0;
for ii = 1:nRP2
	for jj = 1:nMUX
		cnt = cnt+1;
		
		sel		= MX == uMUX(jj) & RP2==uRP2(ii) & abs(AZ)<40 & abs(EL)<40;
		az		= AZ(sel);
		el		= EL(sel);
		chn		= cfg.lookup(sel,1);
		rp2		= RP2(sel);
		mx		= MX(sel);
		
		rpstr = num2str(unique(rp2));
		mxstr = num2str(unique(mx));
		str = ['R' rpstr ', M' mxstr];
		plot(az,el,'ko','MarkerFaceColor',col(cnt,:),'MarkerSize',40);
		text(az,el,str,'HorizontalAlignment','center');
	end
end
xlim([-40 40]);
ylim([-40 40]);
axis square
box off
xlabel('Azimuth (deg)');
ylabel('Elevation (deg)');
title('Zoom');