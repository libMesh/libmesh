clear all
clf
hold on

% These data were gathered using version 1.60 of the CLOC perl script:
% http://sourceforge.net/projects/cloc/files/cloc/v1.60/cloc-1.60.pl/download
%
% Lines of code are counted from include/*/*.h and src/*/*.C.
%
% The basic idea (which should be scripted) is to run:
% git checkout `git rev-list -n 1 --before="$my_date" master`
% cloc.pl src/*/*.C include/*/*.h
%
% The output looks like:
%      799 text files.
%      799 unique files.
%        0 files ignored.
%
% http://cloc.sourceforge.net v 1.60  T=3.43 s (233.1 files/s, 91171.7 lines/s)
% -------------------------------------------------------------------------------
% Language                     files          blank        comment           code
% -------------------------------------------------------------------------------
% C++                            416          43101          34649         135725
% C/C++ Header                   383          20283          42494          36260
% -------------------------------------------------------------------------------
% SUM:                           799          63384          77143         171985
% -------------------------------------------------------------------------------
%
% I only gathered the sum of code lines and the number of files, but
% clearly other data is available.

% Date, Num Files, Num Lines
data = {
% 2012
    '2012-01-04', 694  142345
    '2012-02-04', 697  142585
    '2012-03-04', 703  146127
    '2012-04-04', 706  147191
    '2012-05-04', 708  148202
    '2012-06-04', 705  148334
    '2012-07-04', 713  150066
    '2012-08-04', 727  152269
    '2012-09-04', 725  152381
    '2012-10-04', 1092 155213 % Don't know what happened to the number of files here...
    '2012-11-04', 1094 156082 % We moved from libmesh/src to src around here so maybe that caused it?
    '2012-12-04', 752  156903
% 2013
    '2013-01-04', 754, 158689
    '2013-02-04', 770, 161001
    '2013-03-04', 776, 162189
    '2013-04-04', 783, 162986
    '2013-05-04', 785, 163808
    '2013-06-04', 785, 164022
    '2013-07-04', 789, 163854
    '2013-08-04', 789, 164269
    '2013-09-04', 790, 165129
    '2013-10-04', 790, 165447
    '2013-11-04', 791, 166287
    '2013-12-04', 794, 168772
% 2014
    '2014-01-04', 796, 170174
    '2014-02-04', 796, 170395
    '2014-03-04', 799, 172037
    '2014-04-04', 801, 172262
       };

% length works like you would expect it to for cell arrays.
N=length(data);

% Can't plot string data automatically?
x=linspace(1,N,N);

% Only way to set linestyles with plotyy is after the fact, apparently.
% plotyy returns the axis handle first, followed by the two plot handles.
[haxis, h1, h2] = plotyy(x, [data{:,2}], x, [data{:,3}] ./ 1000);

% Label the axes
ylabel (haxis(1), 'Files');
ylabel (haxis(2), 'Lines of Code (in thousands)');

% Make thick lines - this looks better in PDF
set([h1, h2], 'linewidth', 6);

% Turn on markers
set(h1, 'marker', 'o');
set(h2, 'marker', 's');

% Use bigger markers
set([h1, h2], 'markersize', 12);

% Make dashed line?  Aquaterm doesn't seem to properly display this, but it
% does work in the PDF!
set(h1, 'linestyle', '--');

% Set the xticks and labels so that it shows dates instead of numbers.
xticksat = [1 ceil(N/3) floor(2*N/3) N];
set(haxis, 'xtick', xticksat);

for i=1:length(xticksat)
  xtlabels{i} = data{xticksat(i), 1};
end
set(haxis, 'xticklabel', xtlabels);

% Set some "paper" variables.  None of these seem to have any effect?
% set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]);
% set (gcf, "papersize", [11, 8.5]);
% set (gcf, "paperorientation", 'landscape');

% Make a PDF of this plot.  I had to mess with setting the size
% manually through the print command, otherwise the right-hand axis
% was cut off.  I just picked 1920x1080 at random and it seemed to
% work...
print('-dpdf', '-S1920,1080', 'cloc_libmesh.pdf', '-FHelvetica:20');
