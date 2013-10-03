clear all
clf
hold on

% Hits/month, pages, and gigabytes served.

% To get the Google analytics data:
% .) go to analytics.google.com
% .) click on libmesh
% .) click View Report
% .) Adjust date range to previous month
% .) Record the number of "Pageviews" in the "Hits" column below

%               Hits    , pages,   GB served
cell_data = {
{'Sep', '2013',   9969  ,     0	        , 0.0}
{'Aug', '2013',  10828  ,     0	        , 0.0}
{'Jul', '2013',  14019  ,     0	        , 0.0}
{'Jun', '2013',  13747  ,     0	        , 0.0}
{'May', '2013',  13875  ,     0	        , 0.0}
{'Apr', '2013',  14416  ,     0	        , 0.0}
{'Mar', '2013',  13400  ,     0	        , 0.0}
{'Feb', '2013',  10420  ,     0	        , 0.0}
{'Jan', '2013',  13029  ,     0	        , 0.0}
{'Dec', '2012',  10108  ,     0	        , 0.0} % libmesh switched to github on December 10, 2012
{'Nov', '2012',  11337  ,     0	        , 0.0}
{'Oct', '2012',  13335  ,     0	        , 0.0}
{'Sep', '2012',  13170  ,     0	        , 0.0}
{'Aug', '2012',  13204  ,     0	        , 0.0}
{'Jul', '2012',  12995  ,     0	        , 0.0}
{'Jun', '2012',  12584  ,     0	        , 0.0}
{'May', '2012',  11857  ,     0	        , 0.0}
{'Apr', '2012',  13051  ,     0	        , 0.0}
{'Mar', '2012',  12970  ,     0	        , 0.0}
{'Feb', '2012',  12779  ,     0	        , 0.0}
{'Jan', '2012',  11050  ,     0	        , 0.0}
{'Dec', '2011',  13729  ,     0	        , 0.0}
{'Nov', '2011',  13397  ,     0	        , 0.0}
{'Oct', '2011',  14081  ,     0	        , 0.0}
{'Sep', '2011',  10305  ,     0	        , 0.0}
{'Aug', '2011',  10185  ,     0	        , 0.0} % New "Pageviews" data from google analytics
{'Jul', '2011',      0  ,     0	        , 0.0}
{'Jun', '2011',      0  ,     0	        , 0.0}
{'May', '2011',      0  ,     0	        , 0.0}
{'Apr', '2011',      0  ,     0	        , 0.0}
{'Mar', '2011',      0  ,     0	        , 0.0}
{'Feb', '2011',  34483  , 15002  	, 1.1}
{'Jan', '2011', 133476  , 43549  	, 3.1}
{'Dec', '2010', 105030  , 45570  	, 3.4}
{'Nov', '2010', 115603  , 48695  	, 3.4}
{'Oct', '2010', 98236   , 42091  	, 3.4}
{'Sep', '2010', 81657   , 32305  	, 2.5}
{'Aug', '2010', 99161   , 33199  	, 2.5}
{'Jul', '2010', 87183   , 34026  	, 2.2}
{'Jun', '2010', 128589  , 38019         , 2.5}
{'May', '2010', 117023  , 37107         , 2.5}
{'Apr', '2010', 102493  , 32942         , 2.4}
{'Mar', '2010', 116263  , 42660         , 2.6}
{'Feb', '2010', 208210  , 42680         , 3.0}
{'Jan', '2010', 268856  , 45983         , 3.2}
{'Dec', '2009', 76148   , 37609         , 2.3}
{'Nov', '2009', 88042   , 38869         , 2.5}
{'Oct', '2009', 95727   , 41062         , 2.8}
{'Sep', '2009', 76167   , 33484         , 2.2}
{'Aug', '2009', 64378   , 28542         , 1.9}
{'Jul', '2009', 77828   , 33711         , 2.4}
{'Jun', '2009', 73146   , 31328         , 2.6}
{'May', '2009', 67263   , 26210         , 1.8}
{'Apr', '2009', 70888 	, 26182 	, 1.7 }
{'Mar', '2009', 78200 	, 31139 	, 2.1 }
{'Feb', '2009', 70801 	, 29156 	, 2.1 }
{'Jan', '2009', 73002 	, 29744 	, 2.0 }
{'Dec', '2008', 76248 	, 31690 	, 2.3 }
{'Nov', '2008', 72406 	, 29767 	, 2.3 }
{'Oct', '2008', 64828 	, 27042 	, 2.1 }
{'Sep', '2008', 53057 	, 18635 	, 1.6 }
{'Aug', '2008', 60206 	, 24543 	, 2.0 }
{'Jul', '2008', 60046 	, 22082 	, 2.0 }
{'Jun', '2008', 66055 	, 25418 	, 2.1 }
{'May', '2008', 63123 	, 23277 	, 1.8 }
{'Apr', '2008', 74023 	, 31149 	, 2.0 }
{'Mar', '2008', 69801 	, 29211 	, 1.9 }
{'Feb', '2008', 144881 	, 111751 	, 2.5 }
{'Jan', '2008', 69623 	, 34155 	, 2.0 }
{'Dec', '2007', 80584 	, 42631 	, 2.0 }
{'Nov', '2007', 128307 	, 57851 	, 2.3 }
{'Oct', '2007', 110449 	, 42111 	, 2.4 }
{'Sep', '2007', 143553 	, 39734 	, 2.3 }
{'Aug', '2007', 81192 	, 28187 	, 1.3 }
{'Jul', '2007', 121054 	, 42757 	, 1.8 }
{'Jun', '2007', 115730 	, 38571 	, 1.7 }
{'May', '2007', 182764 	, 129596 	, 1.6 }
{'Apr', '2007', 125578 	, 77539 	, 1.3 }
{'Mar', '2007', 61530 	, 21172 	, 1.2 }
{'Feb', '2007', 52434 	, 17920 	, .9472}
{'Jan', '2007', 53813 	, 19881 	, 1.0  }
{'Dec', '2006', 54265 	, 22688 	, .9162 }
{'Nov', '2006', 55382 	, 19833 	, .9439 }
{'Oct', '2006', 52916 	, 17944 	, .9015 }
{'Sep', '2006', 54919 	, 20591 	, .9056 }
{'Aug', '2006', 55658 	, 11793 	, .6174 }
{'Jul', '2006', 42387 	, 10958 	, .6016 }
{'Jun', '2006', 50495 	, 17128 	, .6881 }
{'May', '2006', 68209 	, 20998 	, .7949 }
{'Apr', '2006', 90314 	, 29635 	, .9762 }
{'Mar', '2006', 47446 	, 14711 	, .6534 }
{'Feb', '2006', 62285 	, 26567 	, .8229 }
{'Jan', '2006', 46723 	, 13487 	, .5662 }
{'Dec', '2005', 90863 	, 40736 	, .9415 }
{'Nov', '2005', 49563 	, 15308 	, .5810 }
{'Oct', '2005', 48336 	, 17976 	, .5749 }
{'Sep', '2005', 57638 	, 14666 	, .4881 }
{'Aug', '2005', 39592 	, 14556 	, .4577 }
{'Jul', '2005', 29944 	, 10006 	, .3363 }
{'Jun', '2005', 40617 	, 16599 	, .5379 }
{'May', '2005', 51288 	, 17629 	, .5689 }
{'Apr', '2005', 29792 	, 8518  	, .2890 }
{'Mar', '2005', 43927 	, 16408 	, .5637 }
{'Feb', '2005', 37733 	, 14077 	, .4373 }
{'Jan', '2005', 3184 	, 3184  	, 0}
{'Dec', '2004', 7483 	, 7483 	        , 0}
{'Nov', '2004', 10261 	, 10261 	, 0}
{'Oct', '2004', 6920 	, 6920 	        , 0}
{'Sep', '2004', 7994 	, 7994 	        , 0}
{'Aug', '2004', 9708 	, 9708 	        , 0}
{'Jul', '2004', 12939 	, 12939 	, 0}
{'Jun', '2004', 12974 	, 12974 	, 0}
{'May', '2004', 11285 	, 11285 	, 0}
{'Apr', '2004', 14995 	, 14995 	, 0}
{'Mar', '2004', 11713 	, 11713 	, 0}
{'Feb', '2004', 11018 	, 11018 	, 0}
{'Jan', '2004', 13599 	, 13599 	, 0}
{'Dec', '2003', 9123 	, 9123 	        , 0}
{'Nov', '2003', 9802 	, 9802 	        , 0}
{'Oct', '2003', 9703 	, 9703 	        , 0}
{'Sep', '2003', 8871 	, 8871 	        , 0}
{'Aug', '2003', 10136 	, 10136 	, 0}
{'Jul', '2003', 6389 	, 6389 	        , 0}
{'Jun', '2003', 6156 	, 6156 	        , 0}
{'May', '2003', 4627 	, 4627 	        , 0}
{'Apr', '2003', 7800 	, 7800 	        , 0}
{'Mar', '2003', 3157 	, 3157 	        , 0}
{'Feb', '2003', 2078 	, 2078 	        , 0}
{'Jan', '2003', 616     , 616           , 0}
};

% Strip out number of hits/month, divided by 1000 for plotting...
N=length(cell_data);
n_hits_month = zeros(N,1);
for i=1:N
  n_hits_month(i) = cell_data{i}{3} / 1000.;
end

% The list is in reverse chronological order,
% so reverse it!
n_hits_month = flipud (n_hits_month);

plot_handle = bar( linspace(1,N,N), n_hits_month, .8 );

% Specify where ticks go in the 'xticksat' array.  Then it
% will set the labels itself.
%           2003 05 07 09 11 13
xticksat = [1    25 49 73 97 121];

% Automatic label setup: put the year
set(gca, 'xtick', xticksat);
for i=1:length(xticksat)
  xtlabels{i} = cell_data{N+1-xticksat(i)}{2};
end
set(gca, 'xticklabel', xtlabels);

% Make sure the ticks stick *out* of the data region
set(gca,'tickdir','out');

% Set a title.  make sure font matches whatever you choose below...
th=title('LibMesh Page Hits/Month (in Thousands)');
set(th, 'fontname', 'Helvetica', 'fontsize', 20);

% Make a PDF of this plot.  This is the only place (I know of)
% where we can control the font size and type used used in the
% X and Y axes...
% orient landscape; % Broken in Octave-3.2.3

% Hard-coded landscape orientation
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]); % [xmin, ymin, xmax, ymax]
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Does this actually work now???  Hurray, it does work as of Octave 3.2.3
set(gca, 'Fontsize', 20);


print('-dpsc', 'libmesh_pagehits.ps', '-FHelvetica:20');
system ('ps2pdf libmesh_pagehits.ps libmesh_pagehits.pdf');
system('rm libmesh_pagehits.ps');



% The old method of gathering this data...

% On the sourceforge site under "Project Web Traffic".  In some
% months, the number of GB served was not recorded: those are marked
% as 0.
%
% This view should give you the last 12 months' data
% https://sourceforge.net/project/stats/detail.php?group_id=71130&ugn=libmesh&mode=12months&type=prweb
%
% For some reason this number jumped to over 260,164 hits for the month of Jan 2010.
% We'll see if it was just a temporary statistics recording error or a change in the
% way they are recording page hits...

% March 2011: Project Web stats are no longer available, please use
% Piwik[1] or Google Analytics[2] to measure your web traffic. Due to site
% security changes, we had to discontinue this part of our stats
% service. Historical stats are still available.
% [1] http://piwik.org/
% [2] http://www.google.com/analytics/
%
% Basically looks like they passed the buck and we just don't have these stats any more...
%

