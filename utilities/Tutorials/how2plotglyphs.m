%% Initialization
close all;
clearvars;


%% Create some data
n	= 10; % number of data points
x	= transpose(1:10);
y	= rand(n,1);


%% Font
% Check which fonts exist.
% Common fonts with useful glyphs:
% - Webdings https://en.wikipedia.org/wiki/Webdings
% - Wingdings https://en.wikipedia.org/wiki/Wingdings
% - Wingdings 2
% - Wingdings 3 
% - ZapfDingbats % https://en.wikipedia.org/wiki/Zapf_Dingbats
% - Zapf Dingbats % https://en.wikipedia.org/wiki/Zapf_Dingbats
listfonts

% We will here plot the loudspeaker glyph from the Webdings Font.
% I picked the Unicode from Adobe Illustrator > Text > Glyphs. Choosing
% Webdings as the font, if you hover above the loudspeaker with 3 waves,
% the unicode U+f055 is shown at the top. This code can be converted to
% uint16 format with the function char.
m			= char(0xF055); % loudspeaker
fontname	= 'Webdings';


%% Plot
figure(1)
clf
h = plot(x,y,'.'); % create an axis to plot in
for ii=1:n
	text(x(ii),y(ii),m,'fontname',fontname,'HorizontalAlignment','center','FontSize',30,'Color','k');
end
delete(h); % delete the dot plot.
xlim([0 11]);
xlabel('x');
ylabel('y');

%% Save as
% Saving as an eps file will remove the fonts from the figure. The way
% to keep the glyphs in a vector-based format, you can export the figure as
% a pdf, through exportgraphics.
exportgraphics(gcf,'vectorglyphfig.pdf')

% The first time I did this, I received an error message about Java Memory.
% This error disappeared when I increased the Java Heap Memory:
% Matlab Home>Preferences>General>Java Heap Memory
