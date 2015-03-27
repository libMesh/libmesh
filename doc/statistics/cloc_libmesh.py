#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num

# git checkout `git rev-list -n 1 --before="$my_date" master`
# cloc.pl src/*/*.C include/*/*.h

data = [
# 2003 - All data from archived svn repo
    # '2003-01-10', 158, 29088, # SVN revision 4 - this is the first revision with trunk/libmesh
    # '2003-01-20', 184, 28937, # SVN revision 11
    # '2003-01-24', 198, 31158, # SVN revision 23
    '2003-02-04', 198, 31344, # SVN revision 47
    '2003-03-04', 243, 36036,
    '2003-04-04', 269, 39946,
    '2003-05-04', 275, 40941,
    '2003-06-04', 310, 44090,
    '2003-07-04', 319, 44445,
    '2003-08-04', 322, 45225,
    '2003-09-04', 325, 46762,
    '2003-10-04', 327, 47151,
    '2003-11-04', 327, 47152, # Up to now, all the include files were in the same directory
    '2003-12-04', 327, 47184,
# 2004 - All data from archived svn repo
    '2004-01-04', 339, 48437,
    '2004-02-04', 343, 50455,
    '2004-03-04', 347, 52198,
    '2004-04-04', 358, 52515,
    '2004-05-04', 358, 52653,
    '2004-06-04', 369, 53953,
    '2004-07-04', 368, 53981,
    '2004-08-04', 371, 54316,
    '2004-09-04', 371, 54510,
    '2004-10-04', 375, 55785,
    '2004-11-04', 375, 55880,
    '2004-12-04', 384, 56612,
# 2005 - All data from archived svn repo
    '2005-01-04', 385, 56674,
    '2005-02-04', 406, 61034,
    '2005-03-04', 406, 62423,
    '2005-04-04', 403, 62595,
    '2005-05-04', 412, 63540,
    '2005-06-04', 416, 69619,
    '2005-07-04', 425, 72092,
    '2005-08-04', 425, 72445,
    '2005-09-04', 429, 74148,
    '2005-10-04', 429, 74263,
    '2005-11-04', 429, 74486,
    '2005-12-04', 429, 74629,
# 2006 - All data from archived svn repo
    '2006-01-04', 429, 74161,
    '2006-02-04', 429, 74165,
    '2006-03-04', 429, 74170,
    '2006-04-04', 429, 74864,
    '2006-05-04', 433, 73847,
    '2006-06-04', 438, 74681,
    '2006-07-04', 454, 76954,
    '2006-08-04', 454, 77464,
    '2006-09-04', 454, 77843,
    '2006-10-04', 454, 78051,
    '2006-11-04', 463, 78683,
    '2006-12-04', 463, 79057,
# 2007 - All data from archived svn repo
    '2007-01-04', 463, 79149,
    '2007-02-04', 475, 79344,
    '2007-03-04', 479, 81416,
    '2007-04-04', 479, 81468,
    '2007-05-04', 481, 84312,
    '2007-06-04', 481, 85565,
    '2007-07-04', 482, 85924,
    '2007-08-04', 485, 86248,
    '2007-09-04', 487, 86481,
    '2007-10-04', 497, 87926,
    '2007-11-04', 502, 89687,
    '2007-12-04', 512, 93523,
# 2008 - All data from archived svn repo
    '2008-01-04', 512, 94263,
    '2008-02-04', 515, 94557,
    '2008-03-04', 526, 98127,
    '2008-04-04', 526, 98256,
    '2008-05-04', 531, 99715,
    '2008-06-04', 531, 99963,
    '2008-07-04', 538, 100839,
    '2008-08-04', 542, 101682,
    '2008-09-04', 548, 102163,
    '2008-10-04', 556, 104185,
    '2008-11-04', 558, 104535,
    '2008-12-04', 565, 106318,
# 2009 - All data from archived svn repo
    '2009-01-04', 565, 106340,
    '2009-02-04', 579, 108431,
    '2009-03-04', 584, 109050,
    '2009-04-04', 584, 109922,
    '2009-05-04', 589, 110821,
    '2009-06-04', 591, 111094,
    '2009-07-04', 591, 111571,
    '2009-08-04', 591, 111555,
    '2009-09-04', 591, 111746,
    '2009-10-04', 591, 111920,
    '2009-11-04', 595, 112993,
    '2009-12-04', 597, 113744,
# 2010 - All data from archived svn repo
    '2010-01-04', 598, 113840,
    '2010-02-04', 600, 114378,
    '2010-03-04', 602, 114981,
    '2010-04-04', 603, 115509,
    '2010-05-04', 603, 115821,
    '2010-06-04', 603, 115875,
    '2010-07-04', 627, 126159,
    '2010-08-04', 627, 126217,
    '2010-09-04', 628, 126078,
    '2010-10-04', 642, 129417,
    '2010-11-04', 643, 130045,
    '2010-12-04', 648, 131363,
# 2011 - All data from archived svn repo
    '2011-01-04', 648, 131644,
    '2011-02-04', 648, 132105,
    '2011-03-04', 658, 132950,
    '2011-04-04', 661, 133643,
    '2011-05-04', 650, 133958,
    '2011-06-04', 662, 134447,
    '2011-07-04', 667, 134938,
    '2011-08-04', 679, 136338,
    '2011-09-04', 684, 138165,
    '2011-10-04', 686, 138627,
    '2011-11-04', 690, 141876,
    '2011-12-04', 690, 142096,
# 2012
    '2012-01-04', 694, 142345,
    '2012-02-04', 697, 142585,
    '2012-03-04', 703, 146127,
    '2012-04-04', 706, 147191,
    '2012-05-04', 708, 148202,
    '2012-06-04', 705, 148334,
    '2012-07-04', 713, 150066,
    '2012-08-04', 727, 152269,
    '2012-09-04', 725, 152381,
    '2012-10-04', 734, 155213, # cloc reports 1092 and 1094 files for Oct/Nov, Don't know what happened...
    '2012-11-04', 743, 156082, # We moved from libmesh/src to src around here so maybe that caused it?
    '2012-12-04', 752, 156903,
# 2013
    '2013-01-04', 754, 158689,
    '2013-02-04', 770, 161001,
    '2013-03-04', 776, 162189,
    '2013-04-04', 783, 162986,
    '2013-05-04', 785, 163808,
    '2013-06-04', 785, 164022,
    '2013-07-04', 789, 163854,
    '2013-08-04', 789, 164269,
    '2013-09-04', 790, 165129,
    '2013-10-04', 790, 165447,
    '2013-11-04', 791, 166287,
    '2013-12-04', 794, 168772,
# 2014
    '2014-01-04', 796, 170174,
    '2014-02-04', 796, 170395,
    '2014-03-04', 799, 172037,
    '2014-04-04', 801, 172262,
    '2014-05-04', 806, 173772,
    '2014-06-04', 807, 171098,
    '2014-07-04', 807, 171220,
    '2014-08-04', 808, 172534,
    '2014-09-04', 808, 173639,
    '2014-10-04', 819, 175750,
    '2014-11-04', 819, 176415,
    '2014-12-04', 819, 176277,
# 2015
    '2015-01-04', 819, 176289,
    '2015-02-04', 824, 176758,
    '2015-03-04', 825, 176958,
       ]

# Extract the dates from the data array
date_strings = data[0::3]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%Y-%m-%d')))

# Extract number of files from data array
n_files = data[1::3]

# Extract number of lines of code from data array
n_lines = data[2::3]

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax1 = fig.add_subplot(111)
ax1.plot(date_nums, n_files, 'bo-')
ax1.set_ylabel('Files (blue circles)')

# Set up x-tick locations
ticks_names = ['2003-03-04', '2007-03-04', '2011-03-04', '2015-01-04']

# Get numerical values for the names
tick_nums = []
for x in ticks_names:
  tick_nums.append(date2num(datetime.strptime(x, '%Y-%m-%d')))

# Set tick labels and positions
ax1.set_xticks(tick_nums)
ax1.set_xticklabels(ticks_names)

# Use the twinx() command to plot more data on the other axis
ax2 = ax1.twinx()
ax2.plot(date_nums, np.divide(n_lines, 1000.), 'gs-')
ax2.set_ylabel('Lines of code in thousands (green squares)')

# Create linear curve fits of the data
files_fit = np.polyfit(date_nums, n_files, 1)
lines_fit = np.polyfit(date_nums, n_lines, 1)

# Convert to files/month
files_per_month = files_fit[0]*(365./12.)
lines_per_month = lines_fit[0]*(365./12.)

# Print curve fit data on the plot , '%.1f'
files_msg = 'Approx. ' + '%.1f' % files_per_month + ' files added/month'
lines_msg = 'Approx. ' + '%.1f' % lines_per_month + ' lines added/month'
ax1.text(date_nums[len(date_nums)/4], 300, files_msg);
ax1.text(date_nums[len(date_nums)/4], 250, lines_msg);

# Save as PDF
plt.savefig('cloc_libmesh.pdf', format='pdf')

# Local Variables:
# python-indent: 2
# End:
