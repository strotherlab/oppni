function print_to_file(figH, outFolder, fname_noext, varargin)

res = '-r350';
format = 'pdf';

set(figH, 'InvertHardCopy', 'off');

% set(figH, 'PaperUnits','centimeters');
% set(figH, 'Units','centimeters');
% 
% pos=get(figH,'Position');
% set(figH, 'PaperSize', [pos(3) pos(4)]);
% % set(figH, 'PaperPositionMode', 'manual');
% set(figH, 'PaperPosition',[0 0 pos(3) pos(4)]);

print(figH, res, ['-d' format], fullfile(outFolder,[ fname_noext '.' format]));

matlab_fig = false;
if nargin > 3
    matlab_fig = varargin{1};
end

if matlab_fig
    saveas(figH, fullfile(outFolder,[ fname_noext '.fig']));
end