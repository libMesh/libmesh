#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num

# SF.net pages and SFLogo Impressions.
# On the site, located under "Sourceforge Traffic".  Number of logos
# (last column) seems to be the most useful one.
#
# This view should give you the last 12 months data
# https://sourceforge.net/project/stats/detail.php?group_id=71130&ugn=libmesh&mode=12months&type=sfweb

# This data has now changed to Google-analytics style...
# After you select the proper date range, scroll down to the bottom
# of the screen and it should show the totals for the two categories,
# which are listed as "SF Logo" and "other"

#               Other   SF Logo
data = [
    'Jan 2003',  681,   479,
    'Feb 2003',  659,  1939,
    'Mar 2003',  488,  1754,
    'Apr 2003',  667,  3202,
    'May 2003',  608,  2552,
    'Jun 2003',  562,  2190,
    'Jul 2003',  997,  3097,
    'Aug 2003',  745,  4708,
    'Sep 2003',  906,  4937,
    'Oct 2003',  892,  6834,
    'Nov 2003', 1257,  8495,
    'Dec 2003', 1147,  6439,
    'Jan 2004',  823,  7791,
    'Feb 2004',  906,  8787,
    'Mar 2004',  979, 11309,
    'Apr 2004',  835,  9393,
    'May 2004',  851,  9796,
    'Jun 2004',  750,  9961,
    'Jul 2004',  618,  6337,
    'Aug 2004',  912,  6647,
    'Sep 2004',  554,  5736,
    'Oct 2004',  524,  6144,
    'Nov 2004',  685,  8122,
    'Dec 2004',  583,  6136,
    'Jan 2005',  215,  2668,
    'Feb 2005',  732,  7618,
    'Mar 2005',  944, 10283,
    'Apr 2005',  837,  9605,
    'May 2005', 1420,  9994,
    'Jun 2005', 1691, 12031,
    'Jul 2005',  849,  6740,
    'Aug 2005', 1068, 11771,
    'Sep 2005', 1119, 11459,
    'Oct 2005',  772,  8614,
    'Nov 2005',  845,  9383,
    'Dec 2005',  814, 10606,
    'Jan 2006', 1004, 11511,
    'Feb 2006',  819, 10693,
    'Mar 2006', 1097, 11925,
    'Apr 2006',  960, 15664,
    'May 2006', 1091, 14194,
    'Jun 2006',  906, 12118,
    'Jul 2006', 1022,  8935,
    'Aug 2006',  914,  9370,
    'Sep 2006', 1087, 11397,
    'Oct 2006', 1311, 11516,
    'Nov 2006', 1182, 10795,
    'Dec 2006',  811,  9418,
    'Jan 2007', 1236, 11522,
    'Feb 2007', 1410, 10669,
    'Mar 2007', 1568, 13141,
    'Apr 2007', 1544, 12285,
    'May 2007', 1362, 14992,
    'Jun 2007', 2229, 17716,
    'Jul 2007', 1822, 15192,
    'Aug 2007', 1446, 12300,
    'Sep 2007', 2045, 19599,
    'Oct 2007', 2680, 14694,
    'Nov 2007', 2344, 15211,
    'Dec 2007', 2235, 10683,
    'Jan 2008', 1582, 11290,
    'Feb 2008', 1712, 12376,
    'Mar 2008', 1908, 13204,
    'Apr 2008', 2308, 13046,
    'May 2008', 2013, 10312,
    'Jun 2008', 2082, 11522,
    'Jul 2008', 1880, 10859,
    'Aug 2008', 2083, 11677,
    'Sep 2008', 1739, 11446,
    'Oct 2008', 2546, 13463,
    'Nov 2008', 2152, 14491,
    'Dec 2008', 2600, 15275,
    'Jan 2009', 1897, 12910,
    'Feb 2009', 1880, 12008,
    'Mar 2009', 6348, 12696,
    'Apr 2009', 1799, 14048,
    'May 2009', 1771, 13122,
    'Jun 2009', 1811, 12114,
    'Jul 2009', 1878, 13600,
    'Aug 2009', 2047, 10828,
    'Sep 2009', 2807, 12914,
    'Oct 2009', 4025, 17326,
    'Nov 2009', 3702, 15648,
    'Dec 2009', 3409, 12510,
    'Jan 2010', 3737, 31211,
    'Feb 2010', 5015, 28772,
    'Mar 2010', 5652, 17882,
    'Apr 2010', 4019, 17495,
    'May 2010', 3336, 18117,
    'Jun 2010', 2174, 21288,
    'Jul 2010',  874, 13900,
    'Aug 2010', 1160, 15153,
    'Sep 2010', 1317, 13836,
    'Oct 2010', 3543, 15279,
    'Nov 2010', 3072, 18663,
    'Dec 2010', 2257, 16381,
    'Jan 2011', 2513, 19798,
    'Feb 2011', 1678, 17870,
    'Mar 2011', 1878, 17759,
    'Apr 2011', 1948, 21264,
    'May 2011', 2696, 15953,
    'Jun 2011', 1514, 18409,
    'Jul 2011', 1422, 13071,
    'Aug 2011',  906,  7857,
    'Sep 2011',  976,  9764,
    'Oct 2011', 1699, 13285,
    'Nov 2011', 1952, 16431,
    'Dec 2011', 2735, 17849,
    'Jan 2012', 1741, 14358,
    'Feb 2012', 1017, 14262,
    'Mar 2012', 1361, 14379,
    'Apr 2012',  967, 15483,
    'May 2012', 2384, 13656,
    'Jun 2012', 1337, 14370,
    'Jul 2012', 2107, 17286,
    'Aug 2012', 8165, 53331,
    'Sep 2012', 2268, 14704,
    'Oct 2012',  738,  7333, # No data recorded from Oct 10 thru 28?
    'Nov 2012', 6104, 39650,
    'Dec 2012', 3439, 24706, # libmesh switched to github Dec 10, 2012
    'Jan 2013', 2552, 31993,
    'Feb 2013', 2107, 24913,
    'Mar 2013', 1376, 23953,
    'Apr 2013', 1582, 19285,
    'May 2013', 1257, 16753,
    'Jun 2013',  482, 14458,
    'Jul 2013',  465, 11325,
    'Aug 2013',  306,  7653,
    'Sep 2013',  731, 11332,
    'Oct 2013',  795, 15619,
    'Nov 2013',  753, 16199,
    'Dec 2013',  593, 11596,
    'Jan 2014',  489, 11195,
    'Feb 2014',  484, 14375,
    'Mar 2014',  363, 13050,
    'Apr 2014',  357, 15700, # As of June 1, 2014 the site above no longer exists...
]

# Extract list of date strings
date_strings = data[0::3]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%b %Y')))

# Strip out number of logos/month for plotting
n_logos_month = data[2::3]

# Scale by 1000
n_logos_month = np.divide(n_logos_month, 1000.)

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Make the bar chart.  One number/month, so width=30
# makes sense.
ax.bar(date_nums, n_logos_month, width=30, color='b')

# Set tick labels at desired locations
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

# Create title
fig.suptitle('SFLogo Pages/Month (in Thousands)')

# Save as PDF
plt.savefig('libmesh_sflogos.pdf')
