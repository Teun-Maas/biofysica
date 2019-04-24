function data_func = getfunc(data_trans)

data_func = data_trans;

%%
trial		= data_func.trial;
ntrial		= numel(trial);

% t			= data_func.time{1};
label		= data_func.label;
labelidx	= vectorstrfind(label,'[O2Hb]'); % find oxy channels. Can we assume deoxy channels are oxy-channel index+1?
nlabel		= numel(labelidx);
label		= strrep(label(labelidx),'O2Hb','functional');
%% Compare oxy and deoxy labels
% labeloxy	= label(idx);
% labelchan	= cellfun(@(x) x(1:end-7),labeloxy,'UniformOutput',false);
% nlabel		= numel(labelchan);
% for ii = 1:nlabel
% 	labelchan{ii}
% end

% sel = t<0;
% close all
kf	= -0.6;
ks	= 0.3;
for ii = 1:ntrial
	trl = trial{ii};
	HbO = trl(labelidx,:);
	HbR = trl(labelidx+1,:);
	
	HbOf	= NaN(size(HbO));
	HbRf		= HbOf;
	for jj = 1:nlabel
		F = functionalvssystemic(HbO(jj,:),HbR(jj,:),kf,ks);
		HbOf(jj,:) = F(1,:);
		HbRf(jj,:) = F(2,:);
	end
	
	trl(labelidx,:) = HbOf;
	trl = trl(labelidx,:);
	% 	trl(labelidx+1,:) = HbRf;
	trial{ii} = trl;
	
	%% Graphics/ test
	% 	subplot(231)
	% 	mu = mean(HbO(:,sel),2);
	% 	x = bsxfun(@minus,HbO,mu)';
	% 	plot(t,x)
	% 	axis square
	%
	%
	% 	subplot(232)
	% 	mu = mean(HbR(:,sel),2);
	% 	y = bsxfun(@minus,HbR,mu)';
	% 	plot(t,y)
	% 	axis square
	%
	% 	subplot(233)
	% 	plot(x,y,'.')
	% 	axis square
	%
	%
	% 	subplot(234)
	% 	mu = mean(HbOf(:,sel),2);
	% 	x = bsxfun(@minus,HbOf,mu)';
	% 	plot(t,x)
	% 	axis square
	%
	%
	% 	subplot(235)
	% 	mu = mean(HbRf(:,sel),2);
	% 	y = bsxfun(@minus,HbRf,mu)';
	% 	plot(t,y)
	% 	axis square
	%
	% 	subplot(236)
	% 	plot(x,y,'.')
	% 	axis square
	
end

data_func.trial = trial;
data_func.label = label;

% save(matfname, 'data_func','-append');