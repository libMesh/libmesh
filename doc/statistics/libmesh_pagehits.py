#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num

# Hits/month, pages, and gigabytes served.

# To get the Google analytics data:
# .) go to analytics.google.com
# .) click on libmesh
# .) click View Report
# .) Adjust date range to previous month
# .) Record the number of "Pageviews" in the "Hits" column below

#                    Hits,  pages, GB served
data = [
#     'Jan 2003',    616,    616, 0
#     'Feb 2003',   2078,   2078, 0,
#     'Mar 2003',   3157,   3157, 0,
#     'Apr 2003',   7800,   7800, 0,
#     'May 2003',   4627,   4627, 0,
#     'Jun 2003',   6156,   6156, 0,
#     'Jul 2003',   6389,   6389, 0,
#     'Aug 2003',  10136,  10136, 0,
#     'Sep 2003',   8871,   8871, 0,
#     'Oct 2003',   9703,   9703, 0,
#     'Nov 2003',   9802,   9802, 0,
#     'Dec 2003',   9123,   9123, 0,
#     'Jan 2004',  13599,  13599, 0,
#     'Feb 2004',  11018,  11018, 0,
#     'Mar 2004',  11713,  11713, 0,
#     'Apr 2004',  14995,  14995, 0,
#     'May 2004',  11285,  11285, 0,
#     'Jun 2004',  12974,  12974, 0,
#     'Jul 2004',  12939,  12939, 0,
#     'Aug 2004',   9708,   9708, 0,
#     'Sep 2004',   7994,   7994, 0,
#     'Oct 2004',   6920,   6920, 0,
#     'Nov 2004',  10261,  10261, 0,
#     'Dec 2004',   7483,   7483, 0,
#     'Jan 2005',   3184,   3184, 0,
#     'Feb 2005',  37733,  14077, .4373,
#     'Mar 2005',  43927,  16408, .5637,
#     'Apr 2005',  29792,   8518, .2890,
#     'May 2005',  51288,  17629, .5689,
#     'Jun 2005',  40617,  16599, .5379,
#     'Jul 2005',  29944,  10006, .3363,
#     'Aug 2005',  39592,  14556, .4577,
#     'Sep 2005',  57638,  14666, .4881,
#     'Oct 2005',  48336,  17976, .5749,
#     'Nov 2005',  49563,  15308, .5810,
#     'Dec 2005',  90863,  40736, .9415,
#     'Jan 2006',  46723,  13487, .5662,
#     'Feb 2006',  62285,  26567, .8229,
#     'Mar 2006',  47446,  14711, .6534,
#     'Apr 2006',  90314,  29635, .9762,
#     'May 2006',  68209,  20998, .7949,
#     'Jun 2006',  50495,  17128, .6881,
#     'Jul 2006',  42387,  10958, .6016,
#     'Aug 2006',  55658,  11793, .6174,
#     'Sep 2006',  54919,  20591, .9056,
#     'Oct 2006',  52916,  17944, .9015,
#     'Nov 2006',  55382,  19833, .9439,
#     'Dec 2006',  54265,  22688, .9162,
#     'Jan 2007',  53813,  19881, 1.0  ,
#     'Feb 2007',  52434,  17920, .9472,
#     'Mar 2007',  61530,  21172, 1.2,
#     'Apr 2007', 125578,  77539, 1.3,
#     'May 2007', 182764, 129596, 1.6,
#     'Jun 2007', 115730,  38571, 1.7,
#     'Jul 2007', 121054,  42757, 1.8,
#     'Aug 2007',  81192,  28187, 1.3,
#     'Sep 2007', 143553,  39734, 2.3,
#     'Oct 2007', 110449,  42111, 2.4,
#     'Nov 2007', 128307,  57851, 2.3,
#     'Dec 2007',  80584,  42631, 2.0,
#     'Jan 2008',  69623,  34155, 2.0,
#     'Feb 2008', 144881, 111751, 2.5,
#     'Mar 2008',  69801,  29211, 1.9,
#     'Apr 2008',  74023,  31149, 2.0,
#     'May 2008',  63123,  23277, 1.8,
#     'Jun 2008',  66055,  25418, 2.1,
#     'Jul 2008',  60046,  22082, 2.0,
#     'Aug 2008',  60206,  24543, 2.0,
#     'Sep 2008',  53057,  18635, 1.6,
#     'Oct 2008',  64828,  27042, 2.1,
#     'Nov 2008',  72406,  29767, 2.3,
#     'Dec 2008',  76248,  31690, 2.3,
#     'Jan 2009',  73002,  29744, 2.0,
#     'Feb 2009',  70801,  29156, 2.1,
#     'Mar 2009',  78200,  31139, 2.1,
#     'Apr 2009',  70888,  26182, 1.7,
#     'May 2009',  67263,  26210, 1.8,
#     'Jun 2009',  73146,  31328, 2.6,
#     'Jul 2009',  77828,  33711, 2.4,
#     'Aug 2009',  64378,  28542, 1.9,
#     'Sep 2009',  76167,  33484, 2.2,
#     'Oct 2009',  95727,  41062, 2.8,
#     'Nov 2009',  88042,  38869, 2.5,
#     'Dec 2009',  76148,  37609, 2.3,
#     'Jan 2010', 268856,  45983, 3.2,
#     'Feb 2010', 208210,  42680, 3.0,
#     'Mar 2010', 116263,  42660, 2.6,
#     'Apr 2010', 102493,  32942, 2.4,
#     'May 2010', 117023,  37107, 2.5,
#     'Jun 2010', 128589,  38019, 2.5,
#     'Jul 2010',  87183,  34026, 2.2,
#     'Aug 2010',  99161,  33199, 2.5,
#     'Sep 2010',  81657,  32305, 2.5,
#     'Oct 2010',  98236,  42091, 3.4,
#     'Nov 2010', 115603,  48695, 3.4,
#     'Dec 2010', 105030,  45570, 3.4,
#     'Jan 2011', 133476,  43549, 3.1,
#     'Feb 2011',  34483,  15002, 1.1,
#     'Mar 2011',      0,      0, 0.0,
#     'Apr 2011',      0,      0, 0.0,
#     'May 2011',      0,      0, 0.0,
#     'Jun 2011',      0,      0, 0.0,
#     'Jul 2011',      0,      0, 0.0,
    'Aug 2011',  10185,      0, 0.0, # New "Pageviews" data from google analytics, does not seem comparable to sf.net pagehits data
    'Sep 2011',  10305,      0, 0.0,
    'Oct 2011',  14081,      0, 0.0,
    'Nov 2011',  13397,      0, 0.0,
    'Dec 2011',  13729,      0, 0.0,
    'Jan 2012',  11050,      0, 0.0,
    'Feb 2012',  12779,      0, 0.0,
    'Mar 2012',  12970,      0, 0.0,
    'Apr 2012',  13051,      0, 0.0,
    'May 2012',  11857,      0, 0.0,
    'Jun 2012',  12584,      0, 0.0,
    'Jul 2012',  12995,      0, 0.0,
    'Aug 2012',  13204,      0, 0.0,
    'Sep 2012',  13170,      0, 0.0,
    'Oct 2012',  13335,      0, 0.0,
    'Nov 2012',  11337,      0, 0.0,
    'Dec 2012',  10108,      0, 0.0, # libmesh switched to github on December 10, 2012
    'Jan 2013',  13029,      0, 0.0,
    'Feb 2013',  10420,      0, 0.0,
    'Mar 2013',  13400,      0, 0.0,
    'Apr 2013',  14416,      0, 0.0,
    'May 2013',  13875,      0, 0.0,
    'Jun 2013',  13747,      0, 0.0,
    'Jul 2013',  14019,      0, 0.0,
    'Aug 2013',  10828,      0, 0.0,
    'Sep 2013',   9969,      0, 0.0,
    'Oct 2013',  13083,      0, 0.0,
    'Nov 2013',  12938,      0, 0.0,
    'Dec 2013',   9079,      0, 0.0,
    'Jan 2014',   9736,      0, 0.0,
    'Feb 2014',  11824,      0, 0.0,
    'Mar 2014',  10861,      0, 0.0,
    'Apr 2014',  12711,      0, 0.0,
    'May 2014',  11177,      0, 0.0,
    'Jun 2014',  10738,      0, 0.0,
    'Jul 2014',  10349,      0, 0.0,
    'Aug 2014',   8877,      0, 0.0,
    'Sep 2014',   9226,      0, 0.0,
    'Oct 2014',   8052,      0, 0.0, # Google analytics number moved over to libmesh.github.io in Oct 2014
    'Nov 2014',   9243,      0, 0.0,
    'Dec 2014',  10714,      0, 0.0,
    'Jan 2015',  11508,      0, 0.0,
    'Feb 2015',  11278,      0, 0.0,
    'Mar 2015',  13305,      0, 0.0,
    'Apr 2015',  12347,      0, 0.0,
    'May 2015',  11368,      0, 0.0,
    'Jun 2015',  11203,      0, 0.0,
    'Jul 2015',  10419,      0, 0.0,
    'Aug 2015',  11282,      0, 0.0,
    'Sep 2015',  13535,      0, 0.0,
    'Oct 2015',  12912,      0, 0.0,
    'Nov 2015',  13894,      0, 0.0,
    'Dec 2015',  11694,      0, 0.0,
    'Jan 2016',  11837,      0, 0.0,
    'Feb 2016',  14102,      0, 0.0,
    'Mar 2016',  13212,      0, 0.0,
    'Apr 2016',  13355,      0, 0.0,
    'May 2016',  12486,      0, 0.0,
    'Jun 2016',  13973,      0, 0.0,
    'Jul 2016',  10688,      0, 0.0,
    'Aug 2016',  10048,      0, 0.0,
    'Sep 2016',  10847,      0, 0.0,
    'Oct 2016',  10984,      0, 0.0,
    'Nov 2016',  12233,      0, 0.0,
    'Dec 2016',  11430,      0, 0.0,
    'Jan 2017',  10327,      0, 0.0,
]

# Extract number of hits/month
n_hits_month = data[1::4]

# Divide by 1000 for plotting...
n_hits_month = np.divide(n_hits_month, 1000.)

# Extract list of date strings
date_strings = data[0::4]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%b %Y')))

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Make the bar chart.  We have one number/month, there are about 30
# days in each month, this defines the bar width...
ax.bar(date_nums, n_hits_month, width=30, color='b')

# Create title
fig.suptitle('LibMesh Page Hits/Month (in Thousands)')

# Set up x-tick locations -- August of each year
ticks_names = ['Aug 2011', 'Aug 2012', 'Aug 2013', 'Aug 2014', 'Aug 2015']

# Get numerical values for the names
tick_nums = []
for x in ticks_names:
  tick_nums.append(date2num(datetime.strptime(x, '%b %Y')))

# Set tick labels and positions
ax.set_xticks(tick_nums)
ax.set_xticklabels(ticks_names)

# Set x limits for the plot
plt.xlim(date_nums[0], date_nums[-1]+30);

# Make x-axis ticks point outward
ax.get_xaxis().set_tick_params(direction='out')

# Save as PDF
plt.savefig('libmesh_pagehits.pdf')

# Local Variables:
# python-indent: 2
# End:
