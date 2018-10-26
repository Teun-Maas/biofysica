function DPF = dpf(lambda,Age,varargin) %#ok<*STOUT>
% D = dpf(wavelength,age)
%
% Molar extinction coefficients for oxygenated and deoxygenated hemoglobin
%
% epsilon:
% - lambda (nm)
% - Hb02  [1/(cm*M)]
% - Hb [1/(cm*M)]
%
% MOLAREXTCOEF('PARAM1',val1) specifies optional
% name/value pairs. Parameters are:
%	'data'		- can be 'default','gratzer','moaveni','takatani' (default
%	= default). The data-set returned are from these people.
%
%	'display'	- display graph. Choices are:
%					0	- no graph (default)
%					>0	- graph
%
% Source:
%
% Scholkmann, F. & Wolf, M. (2013) General equation for the differential
% pathlength factor of the frontal human head depending on wavelength and
% age. J. Biomed. Opt., 18, 105004.
%
% See also: MBLL

% Optional display arguments
disp         = keyval('display',varargin,true);
% dataset      = keyval('data',varargin,'default');
if nargin<1
	lambda = 600:10:800;
end
if nargin<2
	Age = 20:10:50;
end
[lambda,Age] = meshgrid(lambda,Age);


beta(1) = 223.3;
beta(2)  = 0.05624;
beta(3) = 0.8493;
beta(4) = -5.723e-7;
beta(5) = 0.001245;
beta(6)	= -0.9025;

DPF = beta(1)+beta(2)*Age.^beta(3)+beta(4)*lambda.^3+beta(5)*lambda.^2+beta(6)*lambda;

if disp
	graph(DPF,lambda,Age)
end

function graph(DPF,lambda,Age)
close all
figure(1)
clf
subplot(121)
plot(lambda',DPF','-','LineWidth',2)
try
	nicegraph;
catch
end
xlabel('Wavelength (nm)');
ylabel('Differential Pathlength Factor');
legend(num2str(unique(Age)))
subplot(122)
plot(Age,DPF,'-','LineWidth',2)
try
	nicegraph;
catch
end
xlabel('Age (years)');
ylabel('Differential Pathlength Factor');
legend(num2str(unique(lambda)))

% 	hold on
% 	plot(epsilon(:,1),epsilon(:,3),'b-','LineWidth',2)
%
% 	axis square;							% publication-quality
% 	box off;								% publication-quality
% 	set(gca,'TickDir','out','TickLength',[0.005 0.025],...
% 		'XTick',650:50:900,'YTick',500:500:3500,...
% 		'FontSize',15,'FontAngle','Italic','FontName','Helvetica'); % publication-quality
% 	xlabel('Wavelength (nm)','FontSize',15);		% always provide an label on the x-axis
% 	ylabel('Molar extinction coefficient  \epsilon (cm^{-1}M^{-1})','FontSize',15);	% and on the y-axis
% 	axis([600 950 200 4000]);						% set axis limits, this can be maximally between -90 and +90 deg
% 	title('Absorption spectra of hemoglobin'); % Provide all necessary information
% 	legend('oxyhemoglobin','deoxyhemoglobin');

try
	savegraph(mfilename,'eps');
catch
	disp('Why don''t you try out the Biofysica Toolbox @ Gitlab.ru.nl');
end

