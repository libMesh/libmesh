% Number of "papers using libmesh" by year.
%
% Note 1: this does not count citations "only," the authors must have actually
% used libmesh in part of their work.  Therefore, these counts do not include
% things like Wolfgang citing us in his papers to show how Deal.II is
% superior...
%
% Note 2: I typically update this data after regenerating the web page,
% since bibtex2html renumbers the references starting from "1" each year.
%
% Note 3: These citations include anything that is not a dissertation/thesis.
% So, some are conference papers, some are journal articles, etc.
%
% Note 4: The libmesh paper came out in 2006, but there are some citations
% prior to that date, obviously.  These counts include citations of the
% website libmesh.sf.net as well...
%
% Note 5: Preprints are listed as the "current year + 1" and are constantly
% being moved to their respective years after being published.

clear all
clf
hold on

% Global font size for the plot
my_fontsize = 20;

% Use strings to represent the "year"
cell_data = {
    {'2004', 3}
    {'''05', 2}
    {'''06', 11}
    {'''07', 8 }
    {'''08', 22}
    {'''09', 28}
    {'''10', 24}
    {'''11', 35}
    {'''12', 47}
    {'''13', 45}
    {'P', 12} % Preprints
    {'T', 17} % Theses
    };

% Strip numerical data from cell_data
N=length(cell_data);
y=zeros(1,N);
for i=1:N
  y(i) = cell_data{i}{2};
end

% Make bar plot of number of citations
plot_handle = bar( y, .8 );

% Sum up total
total_papers = sum(y);


% Strip x-axis labels from data
xtlabels = cell(1,N);
for i=1:N
  xtlabels(i) = cell_data{i}{1};
end

% Tell plot to use these labels
set(gca, 'xticklabel', xtlabels);

% Where should ticks appear?
xticksat = linspace(1,N,N);
set(gca, 'xtick', xticksat);


% Set the x-min and x-max values
set(gca, 'xlim', [0 N+1]); %data(1,1)-1 data(length(data),1)+1] );

% Make sure the ticks stick *out* of the data region
set(gca,'tickdir','out');

% Set a title.  make sure font matches whatever you choose below...
th=title(['LibMesh Citations, (', num2str(total_papers), ' Total)']);
set(th, 'fontname', 'Helvetica', 'fontsize', my_fontsize);

% Describe the x-axis labels
xh=xlabel('P=Preprints, T=Theses');
set(xh, 'fontname', 'Helvetica', 'fontsize', my_fontsize);

% Hard-coded paper settings... this seemed to be necessary in 3.2.3,
% hopefully they have since fixed it?
set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]); % [xmin, ymin, xmax, ymax]
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% Does this actually work now???  Hurray, it does work as of Octave 3.2.3
set(gca, 'Fontsize', my_fontsize);

print('-dpsc', 'libmesh_citations.ps', '-FHelvetica:20');
system ('ps2pdf libmesh_citations.ps libmesh_citations.pdf');
system('rm libmesh_citations.ps');

