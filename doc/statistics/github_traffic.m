clear all
clf
hold on

% Github has a "traffic" page now, but it doesn't seem like you can
% put in an arbitrary date range?  When I looked at it, it showed the
% numbers of unique visitors and views for the last two weeks only...
% So as long as I check back at least once every two weeks I can keep
% a record of it?  Great...

% https://github.com/libMesh/libmesh/graphs/traffic

% Date, Views, Unique Visitors
data = {
    '2014-Feb-17', 274, 25
    '2014-Feb-18', 145, 30
    '2014-Feb-19', 129, 27
    '2014-Feb-20', 202, 24
    '2014-Feb-21', 240, 22
    '2014-Feb-22', 62,  17
    '2014-Feb-23', 28,  12
    '2014-Feb-24', 217, 19
    '2014-Feb-25', 90,  25
    '2014-Feb-26', 189, 36
    '2014-Feb-27', 112, 26
    '2014-Feb-28', 81,  20
    '2014-Mar-01', 113, 17
    '2014-Mar-02', 53,  16
    '2014-Mar-03', 41,  21
    '2014-Mar-04', 144, 35
    '2014-Mar-05', 51,  20
    '2014-Mar-06', 157, 25
    '2014-Mar-07', 50,  22
    '2014-Mar-08', 50,  11
    '2014-Mar-09', 42,  13
    '2014-Mar-10', 61,  16
    '2014-Mar-11', 27,  16
    '2014-Mar-12', 111, 20
    '2014-Mar-13', 66,  20
    '2014-Mar-14', 223, 25
    '2014-Mar-15', 46,   9
    '2014-Mar-16', 26,  17
    '2014-Mar-17', 80,  29
    '2014-Mar-18', 59,  30
    '2014-Mar-19', 85,  31
    '2014-Mar-20', 122, 18
    '2014-Mar-21', 61,  21
    '2014-Mar-22', 33,  18
    '2014-Mar-23', 64,  14
    '2014-Mar-24', 95,  24
    '2014-Mar-25', 75,  28
    '2014-Mar-26', 49,  18
    '2014-Mar-27', 57,  24
    '2014-Mar-28', 33,  16
    '2014-Mar-29', 41,  16
    '2014-Mar-30', 19,  11
    '2014-Mar-31', 52,  12
    '2014-Apr-01', 120, 21
    '2014-Apr-02', 68,  23
    '2014-Apr-03', 98,  28
    '2014-Apr-04', 77,  21
    '2014-Apr-05', 80,  15
    '2014-Apr-06', 55,  15
    '2014-Apr-07', 71,  31
    '2014-Apr-08', 84,  26
    '2014-Apr-09', 33,  18
    '2014-Apr-10', 32,  16
    '2014-Apr-11', 51,  20
    '2014-Apr-12', 25,  15
    '2014-Apr-13', 49,  20
    };

% length works like you would expect it to for cell arrays.
N=length(data);

% Can't plot string data automatically?
x=linspace(1,N,N);

% Note the extra square brackets!  This seems to be required (at least
% in Octave) to catch the output of data{:,2} as an array.
views = [data{:,2}];
visitors = [data{:,3}];

% Only way to set linestyles with plotyy is after the fact, apparently.
% plotyy returns the axis handle first, followed by the two plot handles.
[haxis, h1, h2] = plotyy(x, views, x, visitors);

% Label the axes
ylabel (haxis(1), 'Daily Page Views');
ylabel (haxis(2), 'Daily Unique Visitors');

% Totaling unique visitors doesn't really mean much...
title (['Total Pageviews: ', num2str(sum(views))]);

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
xticksat = [1 N];
set(haxis, 'xtick', xticksat);

for i=1:length(xticksat)
  xtlabels{i} = data{xticksat(i), 1};
end
set(haxis, 'xticklabel', xtlabels);

% Print to PDF
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]);
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Make a PDF of this plot.
print('-dpdf', 'github_traffic.pdf', '-FHelvetica:20');
