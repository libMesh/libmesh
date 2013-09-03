% This data varies a lot month to month...
% https://sourceforge.net/project/stats/detail.php?group_id=71130&ugn=libmesh&type=svn&mode=12months

clear all
clf
hold on

% After selecting the date range, scroll down to the bottom to see the sums...

% At some point, the names of these data columns changed to:
% .) Read transactions
% .) Write transactions
% .) Write files
% But I think the numbers are still consistent.

% The first month with statistics is October 2007
%               Read    Write   Total files updated
cell_data = {
{'Aug', '2013',  2      ,0      ,0}
{'Jul', '2013',  6      ,0      ,0}
{'Jun', '2013', 34      ,0      ,0}
{'May', '2013', 68      ,0      ,0}
{'Apr', '2013', 80      ,0      ,0}
{'Mar', '2013', 264     ,0      ,0}
{'Feb', '2013', 657     ,0      ,0}
{'Jan', '2013', 10975   ,0      ,0}
{'Dec', '2012', 18856   ,109    ,293} % libmesh switched to github Dec 10, 2012
{'Nov', '2012', 37326   ,637    ,8690}
{'Oct', '2012', 4145    ,185    ,766}
{'Sep', '2012', 125646  ,166    ,3130}
{'Aug', '2012', 9582    ,84     ,279}
{'Jul', '2012', 3722    ,185    ,982}
{'Jun', '2012', 82413   ,121    ,1135}
{'May', '2012', 6936    ,50     ,477}
{'Apr', '2012', 52424   ,146    ,1483}
{'Mar', '2012', 25472   ,338    ,1635}
{'Feb', '2012', 75335   ,180    ,656}
{'Jan', '2012', 17189   ,58     ,4240}
{'Dec', '2011', 3192    ,73     ,1217}
{'Nov', '2011', 39408   ,106    ,346 }
{'Oct', '2011', 58461   ,36     ,110 }
{'Sep', '2011', 6966    ,19     ,42  }
{'Aug', '2011', 71837   ,61     ,230 }
{'Jul', '2011', 46641   ,109    ,2585}
{'Jun', '2011', 3995  	,111  	,459 }
{'May', '2011', 10012  	,71  	,464 }
{'Apr', '2011', 38179  	,128  	,338 }
{'Mar', '2011', 5403  	,105  	,410 }
{'Feb', '2011', 9185  	,30  	,230 }
{'Jan', '2011', 1022  	,28  	,47  }
{'Dec', '2010', 55776  	,27  	,79  }
{'Nov', '2010', 26312  	,42  	,98  }
{'Oct', '2010', 8778  	,64  	,256 }
{'Sep', '2010', 5614  	,58  	,129 }
{'Aug', '2010', 3501  	,64  	,185 }
{'Jul', '2010', 5343    ,15     ,650 }
{'Jun', '2010', 7109    ,50     ,142 }
{'May', '2010', 76137   ,15     ,19  }
{'Apr', '2010', 4048    ,35     ,51  }
{'Mar', '2010', 1004    ,77     ,582 }
{'Feb', '2010', 6846    ,47     ,99  }
{'Jan', '2010', 293     ,37  	,58  }
{'Dec', '2009', 7438  	,27  	,39  }
{'Nov', '2009', 3085 	,46 	,155 }
{'Oct', '2009', 606 	,50 	,274 }
{'Sep', '2009', 3693 	,30 	,51  }
{'Aug', '2009', 3587 	,6 	,13  }
{'Jul', '2009', 374 	,30 	,68  }
{'Jun', '2009', 330 	,22 	,60  }
{'May', '2009', 3560 	,15 	,52  }
{'Apr', '2009', 18652 	,32 	,79  }
{'Mar', '2009', 248 	,22 	,55  }
{'Feb', '2009', 0 	,52 	,178 }
{'Jan', '2009', 0  	,55  	,182 }
{'Dec', '2008', 1  	,29  	,70  }
{'Nov', '2008', 433 	,39 	,116 }
{'Oct', '2008', 305 	,39 	,153 }
{'Sep', '2008', 634 	,72 	,979 }
{'Aug', '2008', 660 	,56 	,103 }
{'Jul', '2008', 144 	,52 	,316 }
{'Jun', '2008', 5656 	,45 	,140 }
{'May', '2008', 147 	,37 	,100 }
{'Apr', '2008', 219 	,53 	,637 }
{'Mar', '2008', 248 	,43 	,102 }
{'Feb', '2008', 5297 	,88 	,240 }
{'Jan', '2008', 110 	,31 	,74  }
{'Dec', '2007', 184  	,44  	,86  }
{'Nov', '2007', 5348 	,258 	,1281}
{'Oct', '2007', 174 	,93 	,633 }
            };

N=length(cell_data);
n_reads   = zeros(N,1);% Probably the most unreliable data...could get reads from anywhere!
n_write   = zeros(N,1);% Misleading since there is no way to tell how many files changed/write
tot_files = zeros(N,1);
for i=1:N
  n_reads(i) = cell_data{i}{3};
  n_write(i) = cell_data{i}{4};
  tot_files(i) = cell_data{i}{5};
end

% The data are recorded in reverse chronological order, so we have to flip them
n_reads   = flipud(n_reads  );
n_write   = flipud(n_write  );
tot_files = flipud(tot_files);

% plot(linspace(1,N,N), n_reads);
% plot(linspace(1,N,N), n_write, 'bo-');
% plot(linspace(1,N,N), tot_files, 'ro-');

plot_handle = bar(linspace(1,N,N), tot_files, .8);

% Specify where ticks go in the 'xticksat' array.  Then it
% will set the labels itself.
%           2008 09 10 11 12 13
xticksat = [4    16 28 40 52 64];

% Set up tick labels on the years.
set(gca, 'xtick', xticksat);
for i=1:length(xticksat)
  % Note: the ordering is reversed by the indexing below...
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
th=title('LibMesh SVN Files Updated/Month');
set(th, 'fontname', 'Helvetica', 'fontsize', 20);

fn='libmesh_svn';
% orient landscape; % Broken in Octave-3.2.3

% Hard-coded landscape orientation
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]); % [xmin, ymin, xmax, ymax]
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Does this actually work now???  Hurray, it does work as of Octave 3.2.3
set(gca, 'Fontsize', 20);

print('-dpsc', [fn,'.ps'], '-FHelvetica:20');
system (['ps2pdf ',fn,'.ps ',fn,'.pdf']);
system(['rm ',fn,'.ps']);

