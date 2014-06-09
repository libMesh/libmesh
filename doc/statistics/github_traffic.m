clear all

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
    '2014-Apr-14', 120, 23
    '2014-Apr-15', 191, 27
    '2014-Apr-16', 219, 24
    '2014-Apr-17', 216, 30
    '2014-Apr-18', 63,  19
    '2014-Apr-19', 36,  11
    '2014-Apr-20', 25,   7
    '2014-Apr-21', 115, 24
    '2014-Apr-22', 128, 31
    '2014-Apr-23', 87,  25
    '2014-Apr-24', 108, 23
    '2014-Apr-25', 111, 20
    '2014-Apr-26', 89,   9
    '2014-Apr-27', 29,  11
    '2014-Apr-28', 177, 28
    '2014-Apr-29', 170, 27
    '2014-Apr-30', 183, 28
    '2014-May-01', 97,  25
    '2014-May-02', 64,  23
    '2014-May-03', 43,  12
    '2014-May-04', 32,  14
    '2014-May-05', 125, 28
    '2014-May-06', 68,  24
    '2014-May-07', 68,  19
    '2014-May-08', 114, 14
    '2014-May-09', 47,  20
    '2014-May-10', 139, 20
    '2014-May-11', 14,   9
    '2014-May-12', 90,  27
    '2014-May-13', 92,  22
    '2014-May-14', 197, 32
    '2014-May-15', 140, 26
    '2014-May-16', 59,  20
    '2014-May-17', 21,  9
    '2014-May-18', 54,  16
    '2014-May-19', 117, 28
    '2014-May-20', 47,  18
    '2014-May-21', 55,  19
    '2014-May-22', 77,  26
    '2014-May-23', 28,  12
    '2014-May-24', 38,  13
    '2014-May-25', 36,  14
    '2014-May-26', 44,  13
    '2014-May-27', 166, 24
    '2014-May-28', 139, 20
    '2014-May-29', 67,  25
    '2014-May-30', 73,  11
    '2014-May-31', 60,  9
    '2014-Jun-01', 22,  11
    '2014-Jun-02', 87,  18
    '2014-Jun-03', 103, 31
    '2014-Jun-04', 105, 27
    '2014-Jun-05', 74,  22
    '2014-Jun-06', 55,  16
    '2014-Jun-07', 53,  15
    '2014-Jun-08', 19,  5
    };

% length works like you would expect it to for cell arrays.
N=length(data);

% A cell array of strings holding just the dates
dates = {data{:,1}};

% Note the extra square brackets!  This seems to be required (at least
% in Octave) to catch the output of data{:,2} as an array.
views = [data{:,2}];
visitors = [data{:,3}];



% 1.) Make daily plot - You can do this, but it gets a bit hard to read

%% % Can't plot string data automatically?
%% x=linspace(1,N,N);
%%
%% % Only way to set linestyles with plotyy is after the fact, apparently.
%% % plotyy returns the axis handle first, followed by the two plot handles.
%% clf;
%% hold on;
%% [haxis, h1, h2] = plotyy(x, views, x, visitors);
%%
%% % Label the axes
%% ylabel (haxis(1), 'Daily Page Views');
%% ylabel (haxis(2), 'Daily Unique Visitors');
%%
%% % Totaling unique visitors doesn't really mean much...
%% title (['Total Pageviews: ', num2str(sum(views)), ', Avg. Daily Unique Visitors: ', num2str(mean(visitors))]);
%%
%% % Make thick lines - this looks better in PDF
%% set([h1, h2], 'linewidth', 6);
%%
%% % Turn on markers
%% set(h1, 'marker', 'o');
%% set(h2, 'marker', 's');
%%
%% % Set marker size
%% set([h1, h2], 'markersize', 6);
%%
%% % Make dashed line?  Aquaterm doesn't seem to properly display this, but it
%% % does work in the PDF!
%% set(h1, 'linestyle', '--');
%%
%% % Set the xticks and labels so that it shows dates instead of numbers.
%% xticksat = [1 length(x)];
%% set(haxis, 'xtick', xticksat);
%%
%% for i=1:length(xticksat)
%%   xtlabels{i} = data{xticksat(i), 1};
%% end
%% set(haxis, 'xticklabel', xtlabels);
%%
%% % Print to PDF
%% set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]);
%% set (gcf, "papersize", [11, 8.5]);
%% set (gcf, "paperorientation", 'landscape');
%%
%% % Make a PDF of this plot.
%% print('-dpdf', 'github_traffic.pdf', '-FHelvetica:20');



% 2.) Make weekly plot - the data is naturally week-periodic
week_indexes = 1:7:N;

% Get total views and average unique viewers for each week
week_views = [];
week_visitors = [];
  for i=1:length(week_indexes)-1
  start = week_indexes(i);
  stop = week_indexes(i+1) - 1;
  week_views(end+1) = sum(views(start:stop));
  week_visitors(end+1) = mean(visitors(start:stop));
end

clf;
hold on;
x=linspace(1, length(week_views), length(week_views));
[haxis, h1, h2] = plotyy(x, week_views, x, week_visitors);

% Label the axes
ylabel (haxis(1), 'Weekly Page Views');
ylabel (haxis(2), 'Avg. Daily Unique Visitors');

% Totaling unique visitors doesn't really mean much...
title (['Total Pageviews: ', num2str(sum(views)), ', Avg. Daily Unique Visitors: ', num2str(mean(visitors))]);

% Make thick lines - this looks better in PDF
set([h1, h2], 'linewidth', 6);

% Turn on markers
set(h1, 'marker', 'o');
set(h2, 'marker', 's');

% Set marker size
set([h1, h2], 'markersize', 6);

% Make dashed line?  Aquaterm doesn't seem to properly display this, but it
% does work in the PDF!
set(h1, 'linestyle', '--');

% Set the xticks and labels so that it shows dates instead of numbers.
xticksat = [1 ceil(length(x)/2) length(x)];
set(haxis, 'xtick', xticksat);

for i=1:length(xticksat)
  xtlabels{i} = ['Week of\n', data{week_indexes(xticksat(i)), 1}];
end
set(haxis, 'xticklabel', xtlabels);

% Print to PDF
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]);
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Make a PDF of this plot.
print('-dpdf', 'weekly_github_traffic.pdf', '-FHelvetica:20');



% 3.) Make monthly plot - not enough data to be interesting yet...

%% month_intervals = {'2014-Feb-17', '2014-Mar-17', '2014-Apr-17'};
%%
%% % Numerical indexes of the month intervals
%% month_indexes = [];
%%
%% for i=1:length(month_intervals)
%%   % Calling 'strcmp' on a cell array of strings returns an array of 0s
%%   % and 1s, then 'find' returns the non-zero index.  The 'end' keyword
%%   % is used to append to the array.
%%   month_indexes(end+1) = find(strcmp(month_intervals{i}, dates));
%% end
%%
%% % Get total views and average unique viewers for each month
%% month_views = [];
%% month_visitors = [];
%%   for i=1:length(month_indexes)-1
%%   start = month_indexes(i);
%%   stop = month_indexes(i+1) - 1;
%%   month_views(end+1) = sum(views(start:stop));
%%   month_visitors(end+1) = mean(visitors(start:stop));
%% end
%%
%% clf;
%% hold on;
%% x=linspace(1, length(month_views), length(month_views));
%% [haxis, h1, h2] = plotyy(x, month_views, x, month_visitors);
