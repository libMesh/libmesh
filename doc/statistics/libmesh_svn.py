#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num

# After selecting the date range, scroll down to the bottom to see the sums...
# At some point, the names of these data columns changed to:
# .) Read transactions
# .) Write transactions
# .) Write files
# But I think the numbers are still consistent.
# If there were no transactions of any kind for a time period (see Sep. 2013)
# sourceforge won't even generate a plot!

# The first month with statistics is October 2007
#               Read    Write   Total files updated
data = [
    'Oct 2007',    174,  93,  633,
    'Nov 2007',   5348, 258, 1281,
    'Dec 2007',    184,  44,   86,
    'Jan 2008',    110,  31,   74,
    'Feb 2008',   5297,  88,  240,
    'Mar 2008',    248,  43,  102,
    'Apr 2008',    219,  53,  637,
    'May 2008',    147,  37,  100,
    'Jun 2008',   5656,  45,  140,
    'Jul 2008',    144,  52,  316,
    'Aug 2008',    660,  56,  103,
    'Sep 2008',    634,  72,  979,
    'Oct 2008',    305,  39,  153,
    'Nov 2008',    433,  39,  116,
    'Dec 2008',      1,  29,   70,
    'Jan 2009',      0,  55,  182,
    'Feb 2009',      0,  52,  178,
    'Mar 2009',    248,  22,   55,
    'Apr 2009',  18652,  32,   79,
    'May 2009',   3560,  15,   52,
    'Jun 2009',    330,  22,   60,
    'Jul 2009',    374,  30,   68,
    'Aug 2009',   3587,   6,   13,
    'Sep 2009',   3693,  30,   51,
    'Oct 2009',    606,  50,  274,
    'Nov 2009',   3085,  46,  155,
    'Dec 2009',   7438,  27,   39,
    'Jan 2010',    293,  37,   58,
    'Feb 2010',   6846,  47,   99,
    'Mar 2010',   1004,  77,  582,
    'Apr 2010',   4048,  35,   51,
    'May 2010',  76137,  15,   19,
    'Jun 2010',   7109,  50,  142,
    'Jul 2010',   5343,  15,  650,
    'Aug 2010',   3501,  64,  185,
    'Sep 2010',   5614,  58,  129,
    'Oct 2010',   8778,  64,  256,
    'Nov 2010',  26312,  42,   98,
    'Dec 2010',  55776,  27,   79,
    'Jan 2011',   1022,  28,   47,
    'Feb 2011',   9185,  30,  230,
    'Mar 2011',   5403, 105,  410,
    'Apr 2011',  38179, 128,  338,
    'May 2011',  10012,  71,  464,
    'Jun 2011',   3995, 111,  459,
    'Jul 2011',  46641, 109, 2585,
    'Aug 2011',  71837,  61,  230,
    'Sep 2011',   6966,  19,   42,
    'Oct 2011',  58461,  36,  110,
    'Nov 2011',  39408, 106,  346,
    'Dec 2011',   3192,  73, 1217,
    'Jan 2012',  17189,  58, 4240,
    'Feb 2012',  75335, 180,  656,
    'Mar 2012',  25472, 338, 1635,
    'Apr 2012',  52424, 146, 1483,
    'May 2012',   6936,  50,  477,
    'Jun 2012',  82413, 121, 1135,
    'Jul 2012',   3722, 185,  982,
    'Aug 2012',   9582,  84,  279,
    'Sep 2012', 125646, 166, 3130,
    'Oct 2012',   4145, 185,  766,
    'Nov 2012',  37326, 637, 8690,
    'Dec 2012',  18856, 109,  293, # libmesh switched to github Dec 10, 2012
    'Jan 2013',  10975,   0,    0,
    'Feb 2013',    657,   0,    0,
    'Mar 2013',    264,   0,    0,
    'Apr 2013',     80,   0,    0,
    'May 2013',     68,   0,    0,
    'Jun 2013',     34,   0,    0,
    'Jul 2013',      6,   0,    0,
    'Aug 2013',      2,   0,    0,
    'Sep 2013',      0,   0,    0,
    'Oct 2013',      0,   0,    0,
    'Nov 2013',      0,   0,    0, # SVN repository deleted from sf.net Nov 11, 2013
    'Dec 2013',      0,   0,    0,
    'Jan 2014',      0,   0,    0,
    'Feb 2014',      0,   0,    0,
    'Mar 2014',      0,   0,    0,
    'Apr 2014',      0,   0,    0, # As of June 1, 2014 the site above no longer exists...
]

# Extract list of date strings
date_strings = data[0::4]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%b %Y')))

# Extract the total number of files updated
tot_files = data[3::4]

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Make the bar chart.
ax.bar(date_nums, tot_files, width=30, color='b')

# Create title
fig.suptitle('LibMesh SVN Files Updated/Month')

# Set tick labels at desired locations
xticklabels = ['Jan\n2008', 'Jan\n2009', 'Jan\n2010', 'Jan\n2011', 'Jan\n2012', 'Jan\n2013']

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
plt.savefig('libmesh_svn.pdf')
