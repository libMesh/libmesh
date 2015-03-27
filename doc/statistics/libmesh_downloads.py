#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num

# As of some time in November 2013, sf switched to google analytics style, and
# I don't think they are reporting MB served anymore.  The easiest way to
# select the previous month is to go to "Specific Date Range" in the drop
# down box...

# Also, when I went back and looked up the October data again, it did not
# agree with what I had previously (155), so it's possible this new data and
# the old data are just not compatible...

# Month, number of d/l, and number of MB served
data = [
    'Jan 2003', 113,  48.2,
    'Feb 2003', 156, 160.5,
    'Mar 2003',  84, 110.6,
    'Apr 2003', 132, 185.2,
    'May 2003',  20,  29.0,
    'Jun 2003',   2,   3.0,
    'Jul 2003', 278, 419.3,
    'Aug 2003', 136, 192.9,
    'Sep 2003', 201, 336.4,
    'Oct 2003', 163, 264.5,
    'Nov 2003', 241, 422.3,
    'Dec 2003', 147, 287.6,
    'Jan 2004', 112, 223.9,
    'Feb 2004', 154, 304.0,
    'Mar 2004', 236, 447.8,
    'Apr 2004', 135, 276.5,
    'May 2004', 172, 351.9,
    'Jun 2004', 223, 446.6,
    'Jul 2004', 129, 262.5,
    'Aug 2004', 199, 398.8,
    'Sep 2004', 109, 226.4,
    'Oct 2004', 146, 302.9,
    'Nov 2004', 170, 357.3,
    'Dec 2004', 121, 247.9,
    'Jan 2005', 118, 246.8,
    'Feb 2005', 120, 249.8,
    'Mar 2005', 179, 378.5,
    'Apr 2005', 133, 280.5,
    'May 2005', 202, 442.6,
    'Jun 2005', 291, 832.0,
    'Jul 2005', 147, 428.8,
    'Aug 2005', 160, 459.8,
    'Sep 2005', 151, 415.8,
    'Oct 2005', 161, 465.6,
    'Nov 2005', 165, 461.7,
    'Dec 2005', 148, 383.5,
    'Jan 2006', 176, 510.4,
    'Feb 2006', 124, 363.4,
    'Mar 2006', 183, 530.2,
    'Apr 2006', 133, 388.4,
    'May 2006', 199, 571.8,
    'Jun 2006', 136, 400.5,
    'Jul 2006', 145, 426.2,
    'Aug 2006', 146, 430.0,
    'Sep 2006', 153, 447.3,
    'Oct 2006', 210, 613.2,
    'Nov 2006', 132, 371.9,
    'Dec 2006', 125, 367.3,
    'Jan 2007', 181, 511.0,
    'Feb 2007', 208, 602.9,
    'Mar 2007', 194, 625.2,
    'Apr 2007', 216, 674.8,
    'May 2007', 201, 650.9,
    'Jun 2007', 291, 923.3,
    'Jul 2007', 236, 816.6,
    'Aug 2007', 201, 700.5,
    'Sep 2007', 224, 744.3,
    'Oct 2007', 328,  1000,
    'Nov 2007', 307,  1000,
    'Dec 2007', 264,  1000,
    'Jan 2008', 178, 700.6,
    'Feb 2008', 177, 684.6,
    'Mar 2008', 233, 911.5,
    'Apr 2008', 231, 895.1,
    'May 2008', 238, 923.3,
    'Jun 2008', 209, 825.4,
    'Jul 2008', 141, 554.5,
    'Aug 2008', 207, 782.3,
    'Sep 2008', 207, 863.7,
    'Oct 2008', 217, 877.1,
    'Nov 2008', 270,  1000,
    'Dec 2008', 222, 839.0,
    'Jan 2009', 153, 600.1,
    'Feb 2009', 265, 785.4,
    'Mar 2009', 226, 743.5,
    'Apr 2009', 145, 617.3,
    'May 2009', 141, 614.6,
    'Jun 2009', 167, 617.8,
    'Jul 2009', 152, 628.7,
    'Aug 2009', 143, 565.2,
    'Sep 2009', 162, 690.6,
    'Oct 2009', 279,  1100,
    'Nov 2009', 285,  1100,
    'Dec 2009', 281,  1000,
    'Jan 2010', 221, 862.7,
    'Feb 2010', 203, 801.8,
    'Mar 2010', 268,  1000,
    'Apr 2010', 175, 767.2,
    'May 2010', 185, 781.1,
    'Jun 2010', 193, 770.1,
    'Jul 2010', 162, 706.9,
    'Aug 2010', 165, 740.3,
    'Sep 2010', 172, 700.3,
    'Oct 2010', 223,  1600,
    'Nov 2010', 266,  1200,
    'Dec 2010', 215,  1000,
    'Jan 2011', 286,  1400,
    'Feb 2011', 172, 943.6,
    'Mar 2011', 216,  1100,
    'Apr 2011', 257,  1300,
    'May 2011', 179, 981.8,
    'Jun 2011', 134, 718.0,
    'Jul 2011', 119, 661.7,
    'Aug 2011',  99,   539,
    'Sep 2011', 162, 827.2,
    'Oct 2011', 203, 843.3,
    'Nov 2011', 210,   0.0,
    'Dec 2011', 216,   0.0,
    'Jan 2012', 169,   0.0,
    'Feb 2012', 199,   0.0,
    'Mar 2012', 230,   0.0,
    'Apr 2012', 243,   0.0,
    'May 2012', 183,   0.0,
    'Jun 2012', 183,   0.0,
    'Jul 2012', 185,   0.0,
    'Aug 2012', 139,   0.0,
    'Sep 2012', 163,   0.0,
    'Oct 2012', 197,   0.0,
    'Nov 2012', 278,   0.0,
    'Dec 2012', 139,   0.0, # This month (around 12/10/2012) libmesh was switched over to github
    'Jan 2013', 108,   0.0,
    'Feb 2013', 250,   0.0,
    'Mar 2013', 259,   0.0,
    'Apr 2013', 282,   0.0,
    'May 2013', 248,   0.0,
    'Jun 2013', 212,   0.0,
    'Jul 2013', 251,   0.0,
    'Aug 2013', 162,   0.0,
    'Sep 2013', 179,   0.0,
    'Oct 2013', 336,   0.0,
    'Nov 2013', 166,   0.0, # All downloadable files removed November 12th
    'Dec 2013',   0,   0.0,
    'Jan 2014',   0,   0.0,
    'Feb 2014',   0,   0.0,
    'Mar 2014',   0,   0.0,
    'Apr 2014',   0,   0.0,
    'May 2014',   0,   0.0,
    'Jun 2014',   0,   0.0,
    'Jul 2014',   0,   0.0,
    'Aug 2014',   0,   0.0,
    'Sep 2014',   0,   0.0,
]

# Extract list of date strings
date_strings = data[0::3]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%b %Y')))

# Strip out number of downloads for plotting
n_downloads = data[1::3]

# Compute the moving average
moving_average = []
running_total = 0
counter = 1
for dl in n_downloads:
  running_total += dl
  moving_average.append(float(running_total) / counter);
  counter += 1

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Plot the moving average
ax.plot(date_nums, moving_average, 'r-', label='All time trailing average');

# Make the bar chart.  We have 1 number/month, so make the width of
# each bar 30 for the approximate number of days in each month.
ax.bar(date_nums, n_downloads, width=30, color='b')

# Add a legend.
plt.legend(loc='upper left')

# Create title
fig.suptitle('LibMesh Downloads/Month')

# Specify tick labels
xticklabels = ['Jan\n2003', 'Jan\n2005', 'Jan\n2007', 'Jan\n2009', 'Jan\n2011', 'Jan\n2013']

# Get numerical values for the tick labels
tick_nums = []
for x in xticklabels:
  tick_nums.append(date2num(datetime.strptime(x, '%b\n%Y')))

ax.set_xticks(tick_nums)
ax.set_xticklabels(xticklabels)

# Make x-axis tick marks point outward
ax.get_xaxis().set_tick_params(direction='out')

# Set the xlimits
plt.xlim(date_nums[0], date_nums[-1]+30);

# Save as PDF
plt.savefig('libmesh_downloads.pdf')

# Local Variables:
# python-indent: 2
# End:
