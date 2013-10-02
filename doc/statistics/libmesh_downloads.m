% This information was copied from sourceforge,
% https://sourceforge.net/project/stats/detail.php?group_id=71130&ugn=libmesh&type=prdownload&mode=alltime&package_id=0
clear all
clf
hold on

% As of some time in November, sf switched to google analytics style, and
% I don't think they are reporting MB served anymore.  The easiest way to
% select the previous month is to go to "Specific Date Range" in the drop
% down box...

% Also, when I went back and looked up the October data again, it did not
% agree with what I had previously (155), so it's possible this new data and
% the old data are just not compatible...

% Month, number of d/l, and number of MB served
cell_data = {
{'Sep', '2013',  	179, 0.0}
{'Aug', '2013',  	162, 0.0}
{'Jul', '2013',  	251, 0.0}
{'Jun', '2013',  	212, 0.0}
{'May', '2013',  	248, 0.0}
{'Apr', '2013',  	282, 0.0}
{'Mar', '2013',  	259, 0.0}
{'Feb', '2013',  	250, 0.0}
{'Jan', '2013',  	108, 0.0}
{'Dec', '2012',  	139, 0.0} % This month (around 12/10/2012) libmesh was switched over to github
{'Nov', '2012',  	278, 0.0}
{'Oct', '2012',  	197, 0.0}
{'Sep', '2012',  	163, 0.0}
{'Aug', '2012',  	139, 0.0}
{'Jul', '2012',  	185, 0.0}
{'Jun', '2012',  	183, 0.0}
{'May', '2012',  	183, 0.0}
{'Apr', '2012',  	243, 0.0}
{'Mar', '2012',  	230, 0.0}
{'Feb', '2012',  	199, 0.0}
{'Jan', '2012',  	169, 0.0}
{'Dec', '2011',  	216, 0.0}
{'Nov', '2011',  	210, 0.0}
{'Oct', '2011',  	203, 843.3}
{'Sep', '2011',  	162, 827.2}
{'Aug', '2011',  	99,  539}
{'Jul', '2011',  	119, 661.7}
{'Jun', '2011',  	134, 718.0}
{'May', '2011',  	179, 981.8}
{'Apr', '2011',  	257, 1300}
{'Mar', '2011',  	216, 1100}
{'Feb', '2011',  	172, 943.6}
{'Jan', '2011',  	286, 1400.0}
{'Dec', '2010',  	215, 1000.0}
{'Nov', '2010',  	266, 1200.0}
{'Oct', '2010',  	223, 1600.0}
{'Sep', '2010',  	172, 700.3}
{'Aug', '2010',  	165, 740.3}
{'Jul', '2010',  	162, 706.9}
{'Jun', '2010',  	193, 770.1}
{'May', '2010',  	185, 781.1}
{'Apr', '2010',  	175, 767.2}
{'Mar', '2010',  	268, 1000.0}
{'Feb', '2010',  	203, 801.8}
{'Jan', '2010',  	221, 862.7}
{'Dec', '2009',  	281, 1000}
{'Nov', '2009',  	285, 1100}
{'Oct', '2009', 	279, 1100}
{'Sep', '2009', 	162, 690.6}
{'Aug', '2009',         143, 565.2 }
{'Jul', '2009',         152, 628.7 }
{'Jun', '2009',         167, 617.8 }
{'May', '2009',         141, 614.6 }
{'Apr', '2009',   	145, 617.3 }
{'Mar', '2009', 	226, 743.5 }
{'Feb', '2009', 	265, 785.4 }
{'Jan', '2009',  	153, 600.1 }
{'Dec', '2008', 	222, 839.0 }
{'Nov', '2008', 	270, 1000.0}
{'Oct', '2008', 	217, 877.1 }
{'Sep', '2008', 	207, 863.7 }
{'Aug', '2008', 	207, 782.3 }
{'Jul', '2008', 	141, 554.5 }
{'Jun', '2008', 	209, 825.4 }
{'May', '2008', 	238, 923.3 }
{'Apr', '2008', 	231, 895.1 }
{'Mar', '2008', 	233, 911.5 }
{'Feb', '2008', 	177, 684.6 }
{'Jan', '2008', 	178, 700.6 }
{'Dec', '2007', 	264, 1000.0}
{'Nov', '2007', 	307, 1000.0}
{'Oct', '2007', 	328, 1000.0}
{'Sep', '2007', 	224, 744.3 }
{'Aug', '2007', 	201, 700.5 }
{'Jul', '2007', 	236, 816.6 }
{'Jun', '2007', 	291, 923.3 }
{'May', '2007', 	201, 650.9 }
{'Apr', '2007', 	216, 674.8 }
{'Mar', '2007', 	194, 625.2 }
{'Feb', '2007', 	208, 602.9 }
{'Jan', '2007', 	181, 511.0 }
{'Dec', '2006', 	125, 367.3 }
{'Nov', '2006', 	132, 371.9 }
{'Oct', '2006', 	210, 613.2 }
{'Sep', '2006', 	153, 447.3 }
{'Aug', '2006', 	146, 430.0 }
{'Jul', '2006', 	145, 426.2 }
{'Jun', '2006', 	136, 400.5 }
{'May', '2006', 	199, 571.8 }
{'Apr', '2006', 	133, 388.4 }
{'Mar', '2006', 	183, 530.2 }
{'Feb', '2006', 	124, 363.4 }
{'Jan', '2006', 	176, 510.4 }
{'Dec', '2005', 	148, 383.5 }
{'Nov', '2005', 	165, 461.7 }
{'Oct', '2005', 	161, 465.6 }
{'Sep', '2005', 	151, 415.8 }
{'Aug', '2005', 	160, 459.8 }
{'Jul', '2005', 	147, 428.8 }
{'Jun', '2005', 	291, 832.0 }
{'May', '2005', 	202, 442.6 }
{'Apr', '2005', 	133, 280.5 }
{'Mar', '2005', 	179, 378.5 }
{'Feb', '2005', 	120, 249.8 }
{'Jan', '2005', 	118, 246.8 }
{'Dec', '2004', 	121, 247.9 }
{'Nov', '2004', 	170, 357.3 }
{'Oct', '2004', 	146, 302.9 }
{'Sep', '2004', 	109, 226.4 }
{'Aug', '2004', 	199, 398.8 }
{'Jul', '2004', 	129, 262.5 }
{'Jun', '2004', 	223, 446.6 }
{'May', '2004', 	172, 351.9 }
{'Apr', '2004', 	135, 276.5 }
{'Mar', '2004', 	236, 447.8 }
{'Feb', '2004', 	154, 304.0 }
{'Jan', '2004', 	112, 223.9 }
{'Dec', '2003', 	147, 287.6 }
{'Nov', '2003', 	241, 422.3 }
{'Oct', '2003', 	163, 264.5 }
{'Sep', '2003', 	201, 336.4 }
{'Aug', '2003', 	136, 192.9 }
{'Jul', '2003', 	278, 419.3 }
{'Jun', '2003', 	2  , 3.0   }
{'May', '2003', 	20 , 29.0  }
{'Apr', '2003', 	132, 185.2 }
{'Mar', '2003', 	84 , 110.6 }
{'Feb', '2003', 	156, 160.5 }
{'Jan', '2003', 	113, 48.2  }
};

% Can't do stuff like...:
%plot(cell_data{:}{1});

% But you can do (both give the same thing...)
% cell_data{1}{1,1}
% cell_data{1}{1}


% Strip out number of downloads for plotting...
N=length(cell_data);
n_downloads = zeros(N,1);
for i=1:N
  n_downloads(i) = cell_data{i}{3};
end

% The list is in reverse chronological order,
% so reverse it!
n_downloads = flipud (n_downloads);

% Plot number of downloads.  The last argument to
% bar is the relative width of each bar.  width=1 means
% there is no whitespace between the bars... Default is
% .8
plot_handle = bar( linspace(1,N,N), n_downloads, .8 );

% Looks like a bar graph with width=1 and with uncolored bars
% plot_handle = stairs( linspace(1,N,N), n_downloads);

% Some interesting properties set by the "bar" command
% facecolor = flat
% facealpha =  1
% You can change the color to black by doing:
% set(plot_handle,'facecolor','k');

% This seems like it should set the transparency but I don't
% think it's supported.
% set(plot_handle,'facealpha',0.5);

% It feels like there is a lot of 'wasted' space in the
% plot when most of the lines are above a certain level.
% We can try to up the minimum y limit... A drawback of
% doing this is that it looks like some data are missing.
%set(gca,'ylim',[100 330]);

% Note, first/last bars need space on either side, so
% don't do this!
% set(gca, 'xlim', [1 N]);

% Compute and plot "all time" moving average
moving_average = zeros(1,N);
running_total = 0;
for i=1:N
  running_total = running_total + n_downloads(i);
  moving_average(i) = running_total / i;
end

ph = plot(linspace(1,N,N), moving_average, 'r-');
set(ph,'linewidth',4); %get(ph) % to see more properties

% Nice, even though this plot was *after* the bar command,
% the lengend goes to the lineplot, not the bar graph!
legend('All time trailing average', 'location', 'northwest');
% Note: 'location', 'outside' is terrible!

% Specify where ticks go in the 'xticksat' array.  Then it
% will set the labels itself.
%           2003 05 07 09 11 13
xticksat = [1    25 49 73 97 121];

% Automatic label setup
set(gca, 'xtick', xticksat);
for i=1:length(xticksat)
  xtlabels{i} = cell_data{N+1-xticksat(i)}{2};
end
set(gca, 'xticklabel', xtlabels);

% Make sure the ticks stick *out* of the data region
set(gca,'tickdir','out');
% set(gca,'ticklength',5); % Unsupported as of 3.0.5, 29th Jan 2010
% set(gca,'linewdith',1.0); % Unsupported  as of 3.0.5, 29th Jan 2010

% When we do this, the extra ticks disappear from the bottom!
% Unfortunately this is not the normal way to read a plot...
% set(gca,'xaxislocation','top');

% Set a title.  make sure font matches whatever you choose below...
th=title('LibMesh Downloads/Month');
set(th, 'fontname', 'Helvetica', 'fontsize', 20);

% Some time around May 2011, the plots all began having rather strange
% xlimits (a big blank space after the the last data).  Setting an explicit
% xlim value seemed to fix that, though.
set(gca, 'xlim', [0 N+1]);

% For some reason orient landscape is broken in octave-3.2.3.  I get
% the warning message:
% warning: print.m - inconsistent papersize and paperorientation properties.
%          papersize = 8.50, 11.00
%          paperorientation = "landscape"
%          the paperorientation property has been ignored

% From the website below, tips for setting paperorientation, papersize, and paperposition manually...
% http://octave.1599824.n4.nabble.com/plot-orientation-feature-plot-filling-page-td1635594.html

% Get current values?  I'm not sure we should bother... just hard-code things because they
% do not seem to be reset even with 'clf'.
%orientation = get (gcf, "paperorientation") % default = 'portrait'
%papersize = get (gcf, "papersize") % default = [8.5 11]
%paperposition = get (gcf, "paperposition") % default = [0.25 2.5 8.0 6.0]

%set (gcf, "paperposition", paperposition([2, 1, 4, 3]));
%set (gcf, "papersize", papersize ([2, 1]));
% setdiff(A,B) returns elements in A that are not in B.
% Do we really need to do this?  What if paperorientation is already correct?
%set (gcf, "paperorientation", setdiff ({"landscape", "portrait"}, {orientation}){1});

% Hard-coded
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]); % [xmin, ymin, xmax, ymax]
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Does this actually work now???  Hurray, it does work as of Octave 3.2.3
set(gca, 'Fontsize', 20);

% Make a PDF of this plot.  This is the only place (I know of)
% where we can control the font size and type used used in the
% X and Y axes... As of Octave 3.2.3, it appears that this method of setting
% the Axes font sizes no longer works either...
print('-dpsc', 'libmesh_downloads.ps', '-FHelvetica:20');
system ('ps2pdf libmesh_downloads.ps libmesh_downloads.pdf');
system('rm libmesh_downloads.ps');
