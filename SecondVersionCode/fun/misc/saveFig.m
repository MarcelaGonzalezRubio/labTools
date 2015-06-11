function saveFig(h,dir,fileName,sizeFlag)
%saveFig saves figure h as .fig and .png for further reference

if nargin<4 || isempty(sizeFlag)
set(h,'Units','Normalized','OuterPosition',[0 0 1 1])
end
fullName=[dir fileName];
if ~exist(dir,'dir')
    mkdir(dir)
end
%print(h, '-painters', '-dpng', '-r900', [fullName '.png']);
saveas(h,[fullName '.fig']) ;
%set(h,'Renderer','opengl');
hgexport(h, [fullName '.png'], hgexport('factorystyle'), 'Format', 'png');

end

