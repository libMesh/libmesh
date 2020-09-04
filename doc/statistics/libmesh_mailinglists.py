#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from operator import add

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num, num2date

# Number of messages to libmesh-devel and libmesh-users over the life
# of the project.  I cut and pasted these from the sf website: they
# are in basically the same format as those pages, which is to say,
# not particularly useful for plotting.

# Administrative interface for both lists.
# https://lists.sourceforge.net/lists/admindb/libmesh-users
# https://lists.sourceforge.net/lists/admindb/libmesh-devel

# You can now see the subscriber counts for both lists in another place:
# https://sourceforge.net/p/libmesh/admin/mailman/

# Month, year, libmesh-devel subscriber count, libmesh-users subscriber count
membership_data = [
    'Jan 2010', 75, 143,
    'Jul 2011', 92, 185,
    'Aug 2011', 91, 184,
    'Sep 2011', 90, 184,
    'Nov 2011', 92, 184,
    'Dec 2011', 93, 185,
    'Jan 2012', 96, 190,
    'Feb 2012', 98, 196,
    'Mar 2012', 99, 199,
    'Apr 2012', 103, 204,
    'May 2012', 105, 210,
    'Jun 2012', 105, 209,
    'Jul 2012', 107, 210,
    'Aug 2012', 109, 213,
    'Sep 2012', 111, 220,
    'Oct 2012', 111, 222,
    'Nov 2012', 112, 225,
    'Dec 2012', 111, 225,
    'Jan 2013', 111, 225,
    'Feb 2013', 112, 228,
    'Mar 2013', 112, 231,
    'Apr 2013', 112, 228,
    'May 2013', 113, 233,
    'Jun 2013', 113, 237,
    'Jul 2013', 114, 240,
    'Aug 2013', 114, 242,
    'Sep 2013', 113, 241,
    'Oct 2013', 112, 241,
    'Nov 2013', 112, 241,
    'Dec 2013', 111, 240,
    'Jan 2014', 112, 244,
    'Feb 2014', 113, 244,
    'Mar 2014', 113, 247,
    'Apr 2014', 113, 248,
    'May 2014', 113, 249,
    'Jun 2014', 113, 247,
    'Jul 2014', 115, 249,
    'Aug 2014', 115, 251,
    'Sep 2014', 117, 254,
    'Oct 2014', 118, 257,
    'Nov 2014', 119, 261,
    'Dec 2014', 120, 262,
    'Jan 2015', 120, 263,
    'Mar 2015', 120, 266, # I missed getting data for Feb 2015
    'Apr 2015', 122, 268,
    'May 2015', 122, 269,
    'Jun 2015', 122, 268,
    'Jul 2015', 122, 271,
    'Aug 2015', 122, 272,
    'Sep 2015', 107, 228,
    'Oct 2015', 107, 230,
    'Nov 2015', 109, 234,
    'Dec 2015', 108, 241,
    'Jan 2016', 108, 239,
    'Feb 2016', 110, 242,
    'Mar 2016', 109, 243,
    'Apr 2016', 110, 242,
    'May 2016', 110, 245,
    'Jun 2016', 110, 244,
    'Jul 2016', 112, 243,
    'Aug 2016', 112, 245,
    'Sep 2016', 115, 246,
    'Oct 2016', 116, 244,
    'Nov 2016', 116, 246,
    'Dec 2016', 116, 246,
    'Jan 2017', 117, 247,
    'Feb 2017', 117, 249,
    'Mar 2017', 117, 250,
    'Apr 2017', 118, 249,
    'May 2017', 118, 249,
    'Jun 2017', 118, 245,
    'Jul 2017', 109, 226,
    'Aug 2017', 109, 226,
    'Sep 2017', 28, 56,   # Sourceforge must have finally dropped the people who didn't reconfirm their mailing list subscriptions!
    'Oct 2017', 30, 63,
    'Nov 2017', 32, 70,
    'Dec 2017', 34, 73,
    'Jan 2018', 35, 77,
    'Feb 2018', 37, 82,
    'Apr 2018', 40, 87,
    'May 2018', 40, 89,
    'Jun 2018', 38, 95,
    'Jul 2018', 39, 99,
    'Aug 2018', 39, 98,
    'Sep 2018', 41, 100,
    'Oct 2018', 41, 103,
    'Nov 2018', 41, 103,
    'Dec 2018', 41, 106,
    'Jan 2019', 41, 104,
    'Feb 2019', 41, 107,
    'Mar 2019', 42, 108,
    'Apr 2019', 42, 107,
    'May 2019', 43, 106,
    'May 2019', 43, 107,
    'Jun 2019', 43, 109,
    'Jul 2019', 42, 107,
    'Aug 2019', 42, 109,
    'Oct 2019', 43, 109,
    'Nov 2019', 43, 111,
    'Dec 2019', 43, 110,
    'Jan 2020', 43, 109,
    'Feb 2020', 44, 114,
    'Mar 2020', 44, 113,
    'Apr 2020', 44, 113,
    'May 2020', 44, 113,
    'Jun 2020', 44, 113,
    'Jul 2020', 44, 113,
    'Aug 2020', 44, 113,
    'Sep 2020', 45, 116,
]

# Strip out the dates from membership_data
date_strings = membership_data[0::3]

# Convert date strings into numbers
date_nums = []
for d in date_strings:
  date_nums.append(date2num(datetime.strptime(d, '%b %Y')))

# Strip out the number of libmesh-devel subscribers from membership_data
devel_count = membership_data[1::3]

# Strip out the number of libmesh-users subscribers from membership_data
users_count = membership_data[2::3]

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# The colors used come from sns.color_palette("muted").as_hex() They
# are the "same basic order of hues as the default matplotlib color
# cycle but more attractive colors."
muted_dark_blue = u'#4878cf'
muted_green = u'#6acc65'
muted_red = u'#d65f5f'
muted_purple = u'#b47cc7'
muted_yellow = u'#c4ad66'
muted_light_blue = u'#77bedb'

# Choose colors from the list above.
primary = muted_dark_blue
secondary = muted_light_blue

# Plot libmesh-users mailing list membership over time
ax.plot(date_nums, users_count, color=primary, marker='s', linestyle='--', label='libmesh-users')

# Plot libmesh-devel mailing list membership over time
ax.plot(date_nums, devel_count, color=secondary, marker='o', linestyle='-', label='libmesh-devel')

# Add a legend to the plot.
plt.legend(loc='upper left')

# Create title
fig.suptitle('LibMesh Mailing List Membership Size')

# Set up xticks and xticklabels
N = len(devel_count)
xtick_indexes = [0, 1, N-1]

xticks = []
xticklabels = []
for index in xtick_indexes:
    xticks.append(date_nums[index])
    xticklabels.append(date_strings[index])

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)

# Save as PDF
plt.savefig('libmesh_mailinglists_membership.pdf')




# libmesh-devel
# https://sourceforge.net/p/libmesh/mailman/libmesh-devel/
#            jan  feb  mar  apr  may  jun  jul  aug  sep  oct  nov  dec
devel_data = [
    '2003',    4,   1,   9,   2,   7,   1,   1,   4,  12,   8,   3,   4,
    '2004',    1,  21,  31,  10,  12,  15,   4,   6,   5,  11,  43,  13,
    '2005',   25,  12,  49,  19, 104,  60,  10,  42,  15,  12,   6,   4,
    '2006',    1,   6,  31,  17,   5,  95,  38,  44,   6,   8,  21,   0,
    '2007',    5,  46,   9,  23,  17,  51,  41,   4,  28,  71, 193,  20,
    '2008',   46,  46,  18,  38,  14, 107,  50, 115,  84,  96, 105,  34,
    '2009',   89,  93, 119,  73,  39,  51,  27,   8,  91,  90,  77,  67,
    '2010',   24,  36,  98,  45,  25,  60,  17,  36,  48,  45,  65,  39,
    '2011',   26,  48, 151, 108,  61, 108,  27,  50,  43,  43,  27,  37,
    '2012',   56, 120,  72,  57,  82,  66,  51,  75, 166, 232, 284, 105, # Dec 10, 2012 libmesh moved to github
    '2013',  168, 151,  30, 145,  26,  53,  76,  33,  23,  72, 125,  38,
    '2014',   47,  62,  27,   8,  12,   2,  22,  22,   0,  17,  20,  12,
    '2015',   25,   2,  16,  13,  21,   5,   1,   8,   9,  30,   8,   0,
    '2016',   16,  31,  43,  18,  21,  11,  17,  26,   4,  16,   5,   6,
    '2017',    1,   2,   5,   4,   1,  11,   5,   0,   3,   1,   7,   0,
    '2018',    8,   8,   1,   0,   5,  11,   0,  51,   3,   0,   0,   0,
    '2019',    2,   0,   3,   7,   2,   0,   6,   0,   0,   4,   0,   0,
    '2020',    0,   0,   0,   0,   1,   0,   0,   0,
]

# libmesh-users starts in Sept 2003!
# https://sourceforge.net/p/libmesh/mailman/libmesh-users/
#            jan  feb  mar  apr  may  jun  jul  aug  sep  oct  nov     dec
users_data = [
    '2003',    0,   0,   0,   0,   0,   0,   0,   0,   2,   2,  27,  31,
    '2004',    6,  15,  33,  10,  46,  11,  21,  15,  13,  23,   1,   8,
    '2005',   27,  57,  86,  23,  37,  34,  24,  17,  50,  24,  10,  60,
    '2006',   47,  46, 127,  19,  26,  62,  47,  51,  61,  42,  50,  33,
    '2007',   60,  55,  77, 102,  82, 102, 169, 117,  80,  37,  51,  43,
    '2008',   71,  94,  98, 125,  54, 119,  60, 111, 118, 125, 119,  94,
    '2009',  109,  38,  93,  88,  29,  57,  53,  48,  68, 151,  23,  35,
    '2010',   84,  60, 184, 112,  60,  90,  23,  70, 119,  27,  47,  54,
    '2011',   22,  19,  92,  93,  35,  91,  32,  61,   7,  69,  81,  23,
    '2012',   64,  95,  35,  36,  63,  98,  70, 171, 149,  64,  67, 126, # Dec 10, 2012 libmesh moved to github
    '2013',  108, 104, 171, 133, 108, 100,  93, 126,  74,  59, 145,  93,
    '2014',   38,  45,  26,  41, 125,  70,  61,  66,  60, 110,  27,  30,
    '2015',   43,  67,  71,  92,  39,  15,  46,  63,  84,  82,  69,  45,
    '2016',   92,  91, 148,  43,  58, 117,  92, 140,  49,  33,  85,  40,
    '2017',   41,  36,  49,  41,  73,  51,  12,  69,  26,  43,  75,  23,
    '2018',   86,  36,  50,  28,  53,  65,  26,  43,  32,  28,  52,  17,
    '2019',   39,  26,  71,  30,  73,  18,   5,  10,   8,  24,  12,  34,
    '2020',   17,  10,   6,   4,  15,   3,   8,  15,
]

# Make plot of monthly data
fig.clf()

# Use a smaller font size on these plots since they are... smaller.
by_month_fontsize = 7

month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Strip out the year strings
year_strings = devel_data[0::13]

for i in range(0, 12):
  # Strip out the devel_data for the current month.  Note that the
  # stride is 13 since the year string also appears in each row.
  devel_data_current_month = devel_data[i+1::13]

  # Strip out the users_data for the current month
  users_data_current_month = users_data[i+1::13]

  # Get the combined number of messages
  combined_data_current_month = np.add(devel_data_current_month, users_data_current_month)

  # Get reference to the ith axes on a 3x4 grid.  Note that the plot
  # numbering starts at 1, simiarly to Matlab.
  ax = fig.add_subplot(3, 4, i+1)

  # The width of the bars to use in the bar chart
  width=.8

  # Make an x-axis to plot against
  N = len(devel_data_current_month)
  x = np.linspace(1, N, N)

  # Plot the summed data
  ax.bar(x, combined_data_current_month, width, color=primary)

  # Plot only the libmesh-devel data
  ax.bar(x, devel_data_current_month, width, color=secondary)

  # Set up the xticks and labels
  xticks = [1, N]
  xticklabels = [year_strings[0], year_strings[N-1]]
  ax.set_xticks(xticks)
  ax.set_xticklabels(xticklabels, fontsize=by_month_fontsize)

  # Set month name as subplot title
  ax.set_title(month_names[i] + ' (max ' + str(max(combined_data_current_month)) + ')', fontsize=by_month_fontsize)

  # Set an empty set of ticks for the y-axis to turn it off.  This
  # is necessary to declutter the figure.
  ax.get_yaxis().set_ticks([])

# We need to leave a bit more room between the subplots to balance the
# font size we can use with the amount of space available on the
# figure.
# wspace = the amount of width reserved for blank space between subplots
# hspace = the amount of height reserved for white space between subplots
plt.subplots_adjust(wspace=.2, hspace=.4)

# Save as PDF
plt.savefig('libmesh_mailinglists_by_month.pdf')



# Make plot of data from all time
fig.clf()

# Strip out all the list entries which are numbers.  Not sure if this
# is the best or most Python way to do this, but it works...
devel_numbers = [x for x in devel_data if (isinstance(x, str)==False)]
users_numbers = [x for x in users_data if (isinstance(x, str)==False)]

# Get the combined number of messages
combined_devel_users_number = np.add(devel_numbers, users_numbers)

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Make an x-axis to plot against
N = len(combined_devel_users_number)
x = np.linspace(1, N, N)

# Plot the combined data
ax.bar(x, combined_devel_users_number, width, color=primary, label='libmesh-users')

# Plot the libmesh-devel data alone
ax.bar(x, devel_numbers, width, color=secondary, label='libmesh-devel')

# Set bi-yearly xticklabels
year_labels = ['2003', '2005', '2007', '2009', '2011', '2013', '2015', '2017', '2019']

# Set up the corresponding tick locations. This starting point was chosen by
# trial and error because it lined up the tick marks fairly well, but I don't
# understand the logic behind it.
xticks = [.55]
for i in range(1, len(year_labels)):
  xticks.append(xticks[i-1] + 24) # 2 years = 24 months

# Center the ticks slightly
xticks = [x+width/2. for x in xticks]
ax.set_xticks(xticks)
ax.set_xticklabels(year_labels)

# Add a legend to the plot.
plt.legend(loc='upper left')

# Set the xlimits
plt.xlim(0, N+2);

# Save as PDF
plt.savefig('libmesh_mailinglists.pdf')

# Local Variables:
# python-indent: 2
# End:
