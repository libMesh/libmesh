#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

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
'Apr', '2014',      0,   0,    0, # As of June 1, 2014 the site above no longer exists...
'Mar', '2014',      0,   0,    0,
'Feb', '2014',      0,   0,    0,
'Jan', '2014',      0,   0,    0,
'Dec', '2013',      0,   0,    0,
'Nov', '2013',      0,   0,    0, # SVN repository deleted from sf.net Nov 11, 2013
'Oct', '2013',      0,   0,    0,
'Sep', '2013',      0,   0,    0,
'Aug', '2013',      2,   0,    0,
'Jul', '2013',      6,   0,    0,
'Jun', '2013',     34,   0,    0,
'May', '2013',     68,   0,    0,
'Apr', '2013',     80,   0,    0,
'Mar', '2013',    264,   0,    0,
'Feb', '2013',    657,   0,    0,
'Jan', '2013',  10975,   0,    0,
'Dec', '2012',  18856, 109,  293, # libmesh switched to github Dec 10, 2012
'Nov', '2012',  37326, 637, 8690,
'Oct', '2012',   4145, 185,  766,
'Sep', '2012', 125646, 166, 3130,
'Aug', '2012',   9582,  84,  279,
'Jul', '2012',   3722, 185,  982,
'Jun', '2012',  82413, 121, 1135,
'May', '2012',   6936,  50,  477,
'Apr', '2012',  52424, 146, 1483,
'Mar', '2012',  25472, 338, 1635,
'Feb', '2012',  75335, 180,  656,
'Jan', '2012',  17189,  58, 4240,
'Dec', '2011',   3192,  73, 1217,
'Nov', '2011',  39408, 106,  346,
'Oct', '2011',  58461,  36,  110,
'Sep', '2011',   6966,  19,   42,
'Aug', '2011',  71837,  61,  230,
'Jul', '2011',  46641, 109, 2585,
'Jun', '2011',   3995, 111,  459,
'May', '2011',  10012,  71,  464,
'Apr', '2011',  38179, 128,  338,
'Mar', '2011',   5403, 105,  410,
'Feb', '2011',   9185,  30,  230,
'Jan', '2011',   1022,  28,   47,
'Dec', '2010',  55776,  27,   79,
'Nov', '2010',  26312,  42,   98,
'Oct', '2010',   8778,  64,  256,
'Sep', '2010',   5614,  58,  129,
'Aug', '2010',   3501,  64,  185,
'Jul', '2010',   5343,  15,  650,
'Jun', '2010',   7109,  50,  142,
'May', '2010',  76137,  15,   19,
'Apr', '2010',   4048,  35,   51,
'Mar', '2010',   1004,  77,  582,
'Feb', '2010',   6846,  47,   99,
'Jan', '2010',    293,  37,   58,
'Dec', '2009',   7438,  27,   39,
'Nov', '2009',   3085,  46,  155,
'Oct', '2009',    606,  50,  274,
'Sep', '2009',   3693,  30,   51,
'Aug', '2009',   3587,   6,   13,
'Jul', '2009',    374,  30,   68,
'Jun', '2009',    330,  22,   60,
'May', '2009',   3560,  15,   52,
'Apr', '2009',  18652,  32,   79,
'Mar', '2009',    248,  22,   55,
'Feb', '2009',      0,  52,  178,
'Jan', '2009',      0,  55,  182,
'Dec', '2008',      1,  29,   70,
'Nov', '2008',    433,  39,  116,
'Oct', '2008',    305,  39,  153,
'Sep', '2008',    634,  72,  979,
'Aug', '2008',    660,  56,  103,
'Jul', '2008',    144,  52,  316,
'Jun', '2008',   5656,  45,  140,
'May', '2008',    147,  37,  100,
'Apr', '2008',    219,  53,  637,
'Mar', '2008',    248,  43,  102,
'Feb', '2008',   5297,  88,  240,
'Jan', '2008',    110,  31,   74,
'Dec', '2007',    184,  44,   86,
'Nov', '2007',   5348, 258, 1281,
'Oct', '2007',    174,  93,  633
            ]

# Extract the total number of files updated
tot_files = data[4::5]

# The list is in reverse chronological order, so reverse it!
tot_files = tot_files[::-1]

# Get a reference to the figure
fig = plt.figure()

# 111 is equivalent to Matlab's subplot(1,1,1) command
ax = fig.add_subplot(111)

# Create an x-axis for plotting
N = len(tot_files);
x = np.linspace(1, N, N)

# Width of the bars
width = 0.8

# Make the bar chart.
ax.bar(x, tot_files, width, color='b')

# Create title
fig.suptitle('LibMesh SVN Files Updated/Month')

# Specify xtick locations
#         2008 09  10  11  12  13
xticks = [4,   16, 28, 40, 52, 64];
xticks = [x+width/2. for x in xticks]
ax.set_xticks(xticks)

# Specify tick labels
xticklabels = ['2008', '2009', '2010', '2011', '2012', '2013']
ax.set_xticklabels(xticklabels)

# Save as PDF
plt.savefig('libmesh_svn.pdf')
