% Number of messages to libmesh-devel and libmesh-users over the life
% of the project.  I cut and pasted these from the sf website: they
% are in basically the same format as those pages, which is to say,
% not particularly useful for plotting.

function libmesh_mailinglists(plot_type)

% plot_type==0: single graph with all data arranged chronologically
% plot_type==1: twelve separate graphs, one for each month.
% plot_type==2: Plot membership data

% clear all % not for functions
clf
hold on

% Approximate subscriber counts:
% https://lists.sourceforge.net/lists/admindb/libmesh-users
% https://lists.sourceforge.net/lists/admindb/libmesh-devel

% Month, year, libmesh-devel subscriber count, libmesh-users subscriber count
membership_data = {
    {'Jan 2010', 75, 143},
    {'Jul 2011', 92, 185},
    {'Aug 2011', 91, 184},
    {'Sep 2011', 90, 184},
    {'Nov 2011', 92, 184},
    {'Dec 2011', 93, 185},
    {'Jan 2012', 96, 190},
    {'Feb 2012', 98, 196},
    {'Mar 2012', 99, 199},
    {'Apr 2012', 103, 204},
    {'May 2012', 105, 210},
    {'Jun 2012', 105, 209},
    {'Jul 2012', 107, 210},
    {'Aug 2012', 109, 213},
    {'Sep 2012', 111, 220},
    {'Oct 2012', 111, 222},
    {'Nov 2012', 112, 225},
    {'Dec 2012', 111, 225},
    {'Jan 2013', 111, 225},
    {'Feb 2013', 112, 228},
    {'Mar 2013', 112, 231}
    {'Apr 2013', 112, 228}
    {'May 2013', 113, 233}
    {'Jun 2013', 113, 237}
    {'Jul 2013', 114, 240}
    {'Aug 2013', 114, 242}
    {'Sep 2013', 113, 241}
    {'Oct 2013', 112, 241}
                  }; % remember to update the indices below!

% The early membership data is spotty, so set indices which are meaningful
% based on the chronology.  Arbitrarily start with index=1.
membership_indices = [1];                                        % 2010
membership_indices = [membership_indices, 12+[7, 8, 9, 11, 12]]; % 2011
membership_indices = [membership_indices, 24+linspace(1,12,12)]; % 2012
membership_indices = [membership_indices, 36+linspace(1,10,10)]; % 2013 <-- Update me

% libmesh-devel
% https://sourceforge.net/mailarchive/forum.php?forum_name=libmesh-devel
devel_data = {
    %        jan        feb     mar     apr     may     jun     jul     aug     sep     oct     nov     dec
    {'2003', 4, 	1, 	9, 	2, 	7, 	1, 	1, 	4, 	12, 	8, 	3, 	4}
    {'2004', 1, 	21, 	31, 	10, 	12, 	15, 	4, 	6, 	5, 	11, 	43, 	13}
    {'2005', 25, 	12, 	49, 	19, 	104, 	60, 	10, 	42, 	15, 	12, 	6, 	4}
    {'2006', 1, 	6, 	31, 	17, 	5, 	95, 	38, 	44, 	6, 	8, 	21, 	0}
    {'2007', 5, 	46, 	9, 	23, 	17, 	51, 	41, 	4, 	28, 	71, 	193, 	20}
    {'2008', 46, 	46, 	18, 	38, 	14, 	107, 	50, 	115, 	84, 	96, 	105, 	34}
    {'2009', 89, 	93, 	119, 	73,     39,     51,     27,     8,      91,     90,     77,     67}
    {'2010', 24,        36,     98,     45,     25,     60,     17,     36,     48,     45,     65,     39}
    {'2011', 26,        48,     151,    108,    61,     108,    27,     50,     43,     43,     27,     37}
    {'2012', 56,        120,    72,     57,     82,     66,     51,     75,    166,    232,    284,    105} % Dec 10, 2012 libmesh moved to github
    {'2013', 168,       151,    30,     145,    26,     53,     76,     33,     23}
    };


% libmesh-users starts in Sept 2003!
% https://sourceforge.net/mailarchive/forum.php?forum_name=libmesh-users
users_data = {
    %           jan     feb     mar     apr     may     jun     jul     aug     sep     oct     nov     dec
    {'2003',    0, 	0, 	0, 	0, 	0, 	0, 	0, 	0, 	2, 	2, 	27, 	31}
    {'2004', 	6, 	15, 	33, 	10, 	46, 	11, 	21, 	15, 	13, 	23, 	1, 	8}
    {'2005', 	27, 	57, 	86, 	23, 	37, 	34, 	24, 	17, 	50, 	24, 	10, 	60}
    {'2006', 	47, 	46, 	127, 	19, 	26, 	62, 	47, 	51, 	61, 	42, 	50, 	33}
    {'2007', 	60, 	55, 	77, 	102, 	82, 	102, 	169, 	117, 	80, 	37, 	51, 	43}
    {'2008', 	71, 	94, 	98, 	125, 	54, 	119, 	60, 	111, 	118, 	125, 	119, 	94}
    {'2009', 	109, 	38, 	93, 	88,     29,     57,     53,     48,     68,     151,    23,     35}
    {'2010',    84,     60,     184,    112,    60,     90,     23,     70,     119,    27,     47,     54}
    {'2011',    22,     19,     92,     93,     35,     91,     32,     61,     7,      69,     81,     23}
    {'2012',    64,     95,     35,     36,     63,     98,     70,     171,    149,    64,     67,    126} % Dec 10, 2012 libmesh moved to github
    {'2013',    108,    104,   171,     133,   108,    100,     93,     126,    74}
    };


if (plot_type==1)
  % Gather data by month
  N=length(devel_data);
  month_names={'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
  bar_x=2003+linspace(0,N-1,N);

  % Set title for all subplots?  This doesn't seem to work...
  % title('libmesh-devel=dk. blue, libmesh-users=cyan');

  for month=1:12
    devel_data_by_month = zeros(1, N);
    devel_plus_users_data_by_month = zeros(1, N);

    for i = 1:N
      % disp(['length(devel_data{i})=', num2str(length(devel_data{i}))]);
      % disp(['month+1=',num2str(month+1)]);
      if (month+1 <= length(devel_data{i}))
        % disp(['Year: ', devel_data{i}{1}, ', Month: ', num2str(month), ', n_msg=', num2str(devel_data{i}{month+1})]);
        devel_data_by_month(i) = devel_data{i}{month+1};
        devel_plus_users_data_by_month(i) = devel_data{i}{month+1} + users_data{i}{month+1};
      end % if
    end % N

    % Go to appropriate subplot and set hold
    subplot(3,4,month);
    hold on;

    % Plot the summed data first
    sum_handle = bar(bar_x, devel_plus_users_data_by_month);
    set(sum_handle,'facecolor','c');

    % Plot the devel-only data second
    dev_handle = bar(bar_x, devel_data_by_month);
    set(dev_handle,'facecolor','b');

    % Find the max summed value to provide scale for each plot
    [max_val, max_yr] = max(devel_plus_users_data_by_month);

    % Label the plot with the Month name using the xlabel.
    % xlabel(month_names{month});

    % Or... set the title as the month name?
    title([month_names{month}, ' (max ', num2str(max_val), ')']);

    % Specify (2) x tick marks
    %xticksat = [bar_x(1) bar_x(N)];

    % Specify (3) x tick marks
    xticksat = [bar_x(1) bar_x(ceil(N/2)) bar_x(N)];

    % Actually Set the tick locations
    set(gca, 'xtick', xticksat);

    % Make sure the ticks stick *out* of the data region
    % This actually isn't great when you are tight on space
    % set(gca,'tickdir','out');

    % Turn off y ticks completely (?)  You will not see any numbers
    % on the axis at all if you do this.
    yticksat = [];
    set(gca, 'ytick', yticksat);

    % Turns on all x tick marks.  Undoes the set(gca, 'xtick', ...) command
    % axis('ticx');


    % Put titles only over certain months?  This doesn't work because
    % it squashes those plots
    % if (month==2)
    %   title('libmesh-devel=dk. blue');
    % end

    % Turn off ytick values?
    % set(gca, 'ytick', []);

    % Use smaller fonts so that they can all fit on the page.
    % This seems to have no effect...
    set(gca, 'Fontsize', 10);

  end % for month

  % A filename to use in the plotting commands
  plot_filename='libmesh_mailinglists_by_month';

  % Set the fontsize to use when calling the 'print' command.
  print_fontsize=10;





elseif (plot_type==0) % single plot over all time

  % Variables used for plotting below...
  plot_filename='libmesh_mailinglists';
  print_fontsize=20;

  % Make arrays with (1) just libmesh-devel and (2) libmesh-devel+libmesh-users.  We will
  % plot the sum as a bar graph behind libmesh-devel using hold on.
  devel_array = [];
  devel_plus_users_array = [];
  ctr = 1;
  for i = 1:length(devel_data)
    for j = 2:length(devel_data{i})
      devel_array(ctr) = devel_data{i}{j};
      devel_plus_users_array(ctr) = devel_data{i}{j} + users_data{i}{j};
      ctr = ctr + 1;
    end
  end

  N=length(devel_array);

  % Sum first
  sum_handle = bar( linspace(1,N,N), devel_plus_users_array );
  set(sum_handle,'facecolor','c');

  % Devel only
  sum_handle = bar( linspace(1,N,N), devel_array );
  set(sum_handle,'facecolor','b');

  % Put up a legend
  legend({'libmesh-users', 'libmesh-devel'}, 'location', 'northwest');


  % Set up x-axis markers.  We'll hard-code this one's years...
  xticksat = [1 25 49 73 97 121];
  xtlabels = {'2003', '2005', '2007', '2009', '2011', '2013'};
  set(gca, 'xtick', xticksat);
  set(gca, 'xticklabel', xtlabels);

  % Make sure the ticks stick *out* of the data region
  set(gca,'tickdir','out');

  % Set a title.  make sure font matches whatever you choose below...
  th=title('LibMesh Mailing List Messages/Month');
  set(th, 'fontname', 'Helvetica', 'fontsize', print_fontsize);

  % Make a PDF of this plot.  This is the only place (I know of)
  % where we can control the font size and type used used in the
  % X and Y axes...
  % orient landscape; % Broken in octave-3.2.3

  % Some time around May 2011, the plots all began having rather strange
  % xlimits (a big blank space after the the last data).  Setting an explicit
  % xlim value seemed to fix that, though.
  set(gca, 'xlim', [0 N+1]);

  % Does this actually work now???  Hurray, it does work as of Octave 3.2.3
  set(gca, 'Fontsize', print_fontsize);




  % Plot list membership data
elseif (plot_type==2)

  % Variables required for plotting below...
  plot_filename='libmesh_mailinglists_membership';
  print_fontsize=20;

  % The amount of data to plot
  N=length(membership_data);

  % Strip out the libmesh devel and users membership numbers
  devel_count = zeros(1, N);
  users_count = zeros(1, N);
  for i=1:N
    devel_count(i) = membership_data{i}{2};
    users_count(i) = membership_data{i}{3};
  end

  ph=plot(membership_indices, users_count, 'cs-');
  set(ph,'linewidth',4);

  hold on;

  ph=plot(membership_indices, devel_count, 'bo-');
  set(ph,'linewidth',4);

  % Put up a legend
  legend({'libmesh-users', 'libmesh-devel'}, 'location', 'northwest');

  % Set a title.
  th=title('LibMesh Mailing List Membership Size');
  set(th, 'fontname', 'Helvetica', 'fontsize', print_fontsize);

  % Set min and max x limits
  set(gca, 'xlim', [0 membership_indices(N)+1]);

  % Set up x-axis markers by choosing where you want the labels to appear in "xtickindices",
  % for example the first, second and Nth data points.
  xtickindices = [1 2 N];

  for i=1:length(xtickindices)
    xticksat(i) = membership_indices( xtickindices(i) );
    xtlabels{i} = membership_data{ xtickindices(i) }{1};
  end

  set(gca, 'xtick', xticksat);
  set(gca, 'xticklabel', xtlabels);


  % You can set axis label font sizes as of Octave 3.2.3
  set(gca, 'Fontsize', print_fontsize);

  % return;
end % if plot_type





% General plotting commands:

% Hard-coded landscape orientation
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]); % [xmin, ymin, xmax, ymax]
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');


print('-dpsc', [plot_filename,'.ps'], ['-FHelvetica:', num2str(print_fontsize)]);
system (['ps2pdf ', plot_filename, '.ps ', plot_filename, '.pdf']);
system(['rm ', plot_filename,'.ps']);

% At some point, Octave started printing directly to PDF so
% print('libmesh_mailinglists.pdf', '-dpdf', '-landscape');

% Unfortunately, the custom fonts don't work quite right yet.  If you
% pass '-FHelvetica' you get an error:
%
% gnuplot> set terminal pdf color"Helvetica"  enhanced;
%                               ^
%         line 0: unexpected text at end of command
% And similarly if you pass '-F:20', GNUplot does not understand it.
%
% According to the GNUplot manual, for the pdf terminal:
% Syntax:
% set terminal pdf {monochrome|color|colour}
%                  {{no}enhanced}
%                  {fname "<font>"} {fsize <fontsize>}
%                  {font "<fontname>{,<fontsize>}"}
%                  {linewidth <lw>} {rounded|butt}
%                  {solid|dashed} {dl <dashlength>}}
%                  {size <XX>{unit},<YY>{unit}}
% There does not seem to be a good way to hack GNUplot to let you do this,
% for example I tried '-F fname Helvetica',
% so we are stuck with the ps terminal still, which does work quite well...


