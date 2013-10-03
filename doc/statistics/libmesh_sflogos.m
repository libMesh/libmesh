clear all
clf
hold on

% SF.net pages and SFLogo Impressions.
% On the site, located under "Sourceforge Traffic".  Number of logos
% (last column) seems to be the most useful one.
%
% This view should give you the last 12 months data
% https://sourceforge.net/project/stats/detail.php?group_id=71130&ugn=libmesh&mode=12months&type=sfweb

% This data has now changed to Google-analytics style...
% After you select the proper date range, scroll down to the bottom
% of the screen and it should show the totals for the two categories,
% which are listed as "SF Logo" and "other"

%               "Other"   "SF Logo"
cell_data = {
{'Sep', '2013',   731   , 11332}
{'Aug', '2013',   306   , 7653}
{'Jul', '2013',   465   , 11325}
{'Jun', '2013',   482   , 14458}
{'May', '2013',  1257   , 16753}
{'Apr', '2013',  1582   , 19285}
{'Mar', '2013',  1376   , 23953}
{'Feb', '2013',  2107   , 24913}
{'Jan', '2013',  2552   , 31993}
{'Dec', '2012',  3439   , 24706} % libmesh switched to github Dec 10, 2012
{'Nov', '2012',  6104   , 39650}
{'Oct', '2012',   738   ,  7333} % No data recorded from Oct 10 thru 28?
{'Sep', '2012',  2268   , 14704}
{'Aug', '2012',  8165   , 53331}
{'Jul', '2012',  2107   , 17286}
{'Jun', '2012',  1337   , 14370}
{'May', '2012',  2384   , 13656}
{'Apr', '2012',   967   , 15483}
{'Mar', '2012',  1361   , 14379}
{'Feb', '2012',  1017   , 14262}
{'Jan', '2012',  1741   , 14358}
{'Dec', '2011',  2735   , 17849}
{'Nov', '2011',  1952	, 16431}
{'Oct', '2011',  1699	, 13285}
{'Sep', '2011',  976	,  9764}
{'Aug', '2011',  906	,  7857}
{'Jul', '2011', 1422	, 13071}
{'Jun', '2011', 1514    , 18409}
{'May', '2011', 2696    , 15953}
{'Apr', '2011', 1948    , 21264}
{'Mar', '2011', 1878    , 17759}
{'Feb', '2011', 1678    , 17870}
{'Jan', '2011', 2513    , 19798}
{'Dec', '2010', 2257    , 16381}
{'Nov', '2010', 3072    , 18663}
{'Oct', '2010', 3543    , 15279}
{'Sep', '2010', 1317    , 13836}
{'Aug', '2010', 1160    , 15153}
{'Jul', '2010', 874     , 13900}
{'Jun', '2010', 2174  	, 21288}
{'May', '2010', 3336  	, 18117}
{'Apr', '2010', 4019  	, 17495}
{'Mar', '2010', 5652  	, 17882}
{'Feb', '2010', 5015   	, 28772}
{'Jan', '2010', 3737   	, 31211}
{'Dec', '2009', 3409   	, 12510}
{'Nov', '2009', 3702   	, 15648}
{'Oct', '2009', 4025  	, 17326}
{'Sep', '2009', 2807  	, 12914}
{'Aug', '2009', 2047    , 10828}
{'Jul', '2009', 1878    , 13600}
{'Jun', '2009', 1811    , 12114}
{'May', '2009', 1771    , 13122}
{'Apr', '2009', 1799  	, 14048}
{'Mar', '2009', 6348 	, 12696}
{'Feb', '2009', 1880 	, 12008}
{'Jan', '2009', 1897 	, 12910}
{'Dec', '2008', 2600 	, 15275}
{'Nov', '2008', 2152 	, 14491}
{'Oct', '2008', 2546 	, 13463}
{'Sep', '2008', 1739 	, 11446}
{'Aug', '2008', 2083 	, 11677}
{'Jul', '2008', 1880 	, 10859}
{'Jun', '2008', 2082 	, 11522}
{'May', '2008', 2013 	, 10312}
{'Apr', '2008', 2308 	, 13046}
{'Mar', '2008', 1908 	, 13204}
{'Feb', '2008', 1712 	, 12376}
{'Jan', '2008', 1582 	, 11290}
{'Dec', '2007', 2235 	, 10683}
{'Nov', '2007', 2344 	, 15211}
{'Oct', '2007', 2680 	, 14694}
{'Sep', '2007', 2045 	, 19599}
{'Aug', '2007', 1446 	, 12300}
{'Jul', '2007', 1822 	, 15192}
{'Jun', '2007', 2229 	, 17716}
{'May', '2007', 1362 	, 14992}
{'Apr', '2007', 1544 	, 12285}
{'Mar', '2007', 1568 	, 13141}
{'Feb', '2007', 1410 	, 10669}
{'Jan', '2007', 1236 	, 11522}
{'Dec', '2006', 811 	, 9418}
{'Nov', '2006', 1182 	, 10795}
{'Oct', '2006', 1311 	, 11516}
{'Sep', '2006', 1087 	, 11397}
{'Aug', '2006', 914 	, 9370}
{'Jul', '2006', 1022 	, 8935}
{'Jun', '2006', 906 	, 12118}
{'May', '2006', 1091 	, 14194}
{'Apr', '2006', 960 	, 15664}
{'Mar', '2006', 1097 	, 11925}
{'Feb', '2006', 819 	, 10693}
{'Jan', '2006', 1004 	, 11511}
{'Dec', '2005', 814 	, 10606}
{'Nov', '2005', 845 	, 9383}
{'Oct', '2005', 772 	, 8614}
{'Sep', '2005', 1119 	, 11459}
{'Aug', '2005', 1068 	, 11771}
{'Jul', '2005', 849 	, 6740}
{'Jun', '2005', 1691 	, 12031}
{'May', '2005', 1420 	, 9994}
{'Apr', '2005', 837 	, 9605}
{'Mar', '2005', 944 	, 10283}
{'Feb', '2005', 732 	, 7618}
{'Jan', '2005', 215 	, 2668}
{'Dec', '2004', 583 	, 6136}
{'Nov', '2004', 685 	, 8122}
{'Oct', '2004', 524 	, 6144}
{'Sep', '2004', 554 	, 5736}
{'Aug', '2004', 912 	, 6647}
{'Jul', '2004', 618 	, 6337}
{'Jun', '2004', 750 	, 9961}
{'May', '2004', 851 	, 9796}
{'Apr', '2004', 835 	, 9393}
{'Mar', '2004', 979 	, 11309}
{'Feb', '2004', 906 	, 8787}
{'Jan', '2004', 823 	, 7791}
{'Dec', '2003', 1147 	, 6439}
{'Nov', '2003', 1257 	, 8495}
{'Oct', '2003', 892 	, 6834}
{'Sep', '2003', 906 	, 4937}
{'Aug', '2003', 745 	, 4708}
{'Jul', '2003', 997 	, 3097}
{'Jun', '2003', 562 	, 2190}
{'May', '2003', 608 	, 2552}
{'Apr', '2003', 667 	, 3202}
{'Mar', '2003', 488 	, 1754}
{'Feb', '2003', 659 	, 1939}
{'Jan', '2003', 681 	, 479}
};



% Strip out number of hits/month, divided by 1000for plotting...
N=length(cell_data);
n_logos_month = zeros(N,1);
for i=1:N
  n_logos_month(i) = cell_data{i}{4} / 1000.;
end

% The list is in reverse chronological order,
% so reverse it!
n_logos_month = flipud (n_logos_month);

plot_handle = bar( linspace(1,N,N), n_logos_month, .8 );

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

% Some time around May 2011, the plots all began having rather strange
% xlimits (a big blank space after the the last data).  Setting an explicit
% xlim value seemed to fix that, though.
set(gca, 'xlim', [0 N+1]);

% Set a title.  make sure font matches whatever you choose below...
th=title('SFLogo Pages/Month (in Thousands)');
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

print('-dpsc', 'libmesh_sflogos.ps', '-FHelvetica:20');
system ('ps2pdf libmesh_sflogos.ps libmesh_sflogos.pdf');
system('rm libmesh_sflogos.ps');
