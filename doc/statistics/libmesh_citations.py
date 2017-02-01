#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Number of "papers using libmesh" by year.
#
# Note 1: this does not count citations "only," the authors must have actually
# used libmesh in part of their work.  Therefore, these counts do not include
# things like Wolfgang citing us in his papers to show how Deal.II is
# superior...
#
# Note 2: I typically update this data after regenerating the web page,
# since bibtex2html renumbers the references starting from "1" each year.
#
# Note 3: These citations include anything that is not a dissertation/thesis.
# So, some are conference papers, some are journal articles, etc.
#
# Note 4: The libmesh paper came out in 2006, but there are some citations
# prior to that date, obviously.  These counts include citations of the
# website libmesh.sf.net as well...
#
# Note 5: Preprints are listed as the "current year + 1" and are constantly
# being moved to their respective years after being published.

data = [
    '2004',  5,
    '\'05',  2,
    '\'06', 13,
    '\'07',  8,
    '\'08', 24,
    '\'09', 30,
    '\'10', 25,
    '\'11', 37,
    '\'12', 50,
    '\'13', 79,
    '\'14', 60,
    '\'15', 78,
    '\'16', 81,
    '\'17', 13,
    'P',     6, # Preprints
    'T',    57  # Theses
    ]

# Extract the x-axis labels from the data array
xlabels = data[0::2]

# Extract the publication counts from the data array
n_papers = data[1::2]

# The number of data points
N = len(xlabels);

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Create an x-axis for plotting
x = np.linspace(1, N, N)

# Width of the bars
width = 0.8

# Make the bar chart.  Plot years in blue, preprints and theses in green.
ax.bar(x[0:N-2], n_papers[0:N-2], width, color='b')
ax.bar(x[N-2:N], n_papers[N-2:N], width, color='g')

# Label the x-axis
plt.xlabel('P=Preprints, T=Theses')

# Set up the xtick locations and labels.  Note that you have to offset
# the position of the ticks by width/2, where width is the width of
# the bars.
ax.set_xticks(np.linspace(1,N,N) + width/2)
ax.set_xticklabels(xlabels)

# Create a title string
title_string = 'Papers by People Using LibMesh, (' + str(sum(n_papers)) + ' Total)'
fig.suptitle(title_string)

# Save as PDF
plt.savefig('libmesh_citations.pdf')

# Local Variables:
# python-indent: 2
# End:
