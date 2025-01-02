#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import datetime
from matplotlib.dates import date2num, num2date
import argparse

rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 13
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True

# Small font size for the legend
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('xx-small')

# These colors used come from:
#
# pip install seaborn
# import seaborn as sns
# sns.color_palette("muted").as_hex()
#
# They are the "same basic order of hues as the default matplotlib
# color cycle but more attractive colors."
muted_dark_blue = u'#4878cf'
muted_green = u'#6acc65'
muted_red = u'#d65f5f'
muted_purple = u'#b47cc7'
muted_yellow = u'#c4ad66'
muted_light_blue = u'#77bedb'
muted_orange = u'#ee854a'
muted_brown = u'#8c613c'
muted_pink = u'#dc7ec0'
muted_grey = u'#797979'
muted_gold = u'#d5bb67'

fig = plt.figure()
ax1 = fig.add_subplot(111)

"""
Plot a running total of the number of publications over the course of a year,
plus some overlap on each end, since publications from year N often come out
in year N-1 or N+1.

Remarks:
* At some point in early 2020 I moved theses and dissertations into
  their respective year files, so there is a relatively large jump
  around that time.
* It is possible for the sum to decrease if a preprint is released
  in year N but then doesn't make it to "print" until year N+1.
* To check out a revision on the first of a month:
  git checkout `git rev-list -n 1 --before="2018-01-01" master`
"""

# 2017 monthly publication totals (first of each month)
# These numbers don't include the number of theses/dissertations.
pub_2017 = [
'2016-10-01', 0,
'2016-11-01', 3,
'2016-12-01', 5,
'2017-01-01', 7, # new year
'2017-02-01', 13,
'2017-03-01', 19,
'2017-04-01', 24,
'2017-05-01', 31,
'2017-06-01', 41,
'2017-07-01', 47,
'2017-08-01', 64,
'2017-09-01', 70,
'2017-10-01', 89,
'2017-11-01', 104,
'2017-12-01', 117,
'2018-01-01', 118, # new year
'2018-02-01', 119,
'2018-03-01', 119,
'2018-04-01', 118,
]

# 2018 monthly publication totals (first of each month)
pub_2018 = [
'2017-10-01', 0,
'2017-11-01', 4,
'2017-12-01', 10,
'2018-01-01', 14, # new year
'2018-02-01', 20,
'2018-03-01', 27,
'2018-04-01', 35,
'2018-05-01', 42,
'2018-06-01', 49,
'2018-07-01', 57,
'2018-08-01', 74,
'2018-09-01', 89,
'2018-10-01', 103,
'2018-11-01', 118,
'2018-12-01', 136,
'2019-01-01', 141, # new year
'2019-02-01', 146,
'2019-03-01', 150,
'2019-04-01', 149,
]

# 2019 monthly publication totals (first of each month)
pub_2019 = [
  '2018-10-01', 0,
  '2018-11-01', 3,
  '2018-12-01', 9,
  '2019-01-01', 11, # new year
  '2019-02-01', 23,
  '2019-03-01', 31,
  '2019-04-01', 44,
  '2019-05-01', 46,
  '2019-06-01', 61,
  '2019-07-01', 72,
  '2019-08-01', 77,
  '2019-09-01', 85,
  '2019-10-01', 91,
  '2019-11-01', 98,
  '2019-12-01', 115,
  '2020-01-01', 119, # new year
  '2020-02-01', 123,
  '2020-03-01', 133,
  '2020-04-01', 132,
]

# 2020 monthly publication totals (first of each month)
pub_2020 = [
  '2019-08-01', 0,
  '2019-09-01', 1,
  '2019-10-01', 3,
  '2019-11-01', 6,
  '2019-12-01', 11,
  '2020-01-01', 17, # new year
  '2020-02-01', 34,
  '2020-03-01', 44,
  '2020-04-01', 57,
  '2020-05-01', 75,
  '2020-06-01', 96,
  '2020-07-01', 111,
  '2020-08-01', 122,
  '2020-09-01', 139,
  '2020-10-01', 144,
  '2020-11-01', 158,
  '2020-12-01', 178,
  '2021-01-01', 187, # new year
  '2021-02-01', 189,
  '2021-03-01', 191,
  '2021-04-01', 196,
]

# 2021 monthly publication totals (first of each month)
pub_2021 = [
  '2020-09-01', 0,
  '2020-10-01', 1,
  '2020-11-01', 4,
  '2020-12-01', 8,
  '2021-01-01', 11, # new year
  '2021-02-01', 20,
  '2021-03-01', 30,
  '2021-04-01', 44,
  '2021-05-01', 53,
  '2021-06-01', 74,
  '2021-07-01', 82,
  '2021-08-01', 113,
  '2021-09-01', 121,
  '2021-10-01', 144,
  '2021-11-01', 162,
  '2021-12-01', 173,
  '2022-01-01', 203, # new year
  '2022-02-01', 203,
  '2022-03-01', 204,
  '2022-04-01', 202,
]

# 2022 monthly publication totals (first of each month)
pub_2022 = [
  '2021-09-01', 0,
  '2021-10-01', 2,
  '2021-11-01', 6,
  '2021-12-01', 9,
  '2022-01-01', 13, # new year
  '2022-02-01', 28,
  '2022-03-01', 47,
  '2022-04-01', 66,
  '2022-05-01', 78,
  '2022-06-01', 101,
  '2022-07-01', 117,
  '2022-08-01', 136,
  '2022-09-01', 151,
  '2022-10-01', 172,
  '2022-11-01', 195,
  '2022-12-01', 211,
  '2023-01-01', 216, # new year
  '2023-02-01', 220,
  '2023-03-01', 222,
  '2023-04-01', 221,
]

# 2023 monthly publication totals (first of each month)
pub_2023 = [
  '2022-06-01', 0,
  '2022-07-01', 1,
  '2022-08-01', 1,
  '2022-09-01', 1,
  '2022-10-01', 1,
  '2022-11-01', 2,
  '2022-12-01', 6,
  '2023-01-01', 14, # new year
  '2023-02-01', 27,
  '2023-03-01', 44,
  '2023-04-01', 65,
  '2023-05-01', 71,
  '2023-06-01', 93,
  '2023-07-01', 109,
  '2023-08-01', 137,
  '2023-09-01', 161,
  '2023-10-01', 189,
  '2023-11-01', 206,
  '2023-12-01', 226,
  '2024-01-01', 232, # new year
  '2024-02-01', 235,
  '2024-03-01', 234,
  '2024-04-01', 237,
]

# 2024 monthly publication totals (first of each month)
pub_2024 = [
  '2023-10-01', 0,
  '2023-11-01', 1,
  '2023-12-01', 3,
  '2024-01-01', 14, # new year
  '2024-02-01', 35,
  '2024-03-01', 50,
  '2024-04-01', 67,
  '2024-05-01', 85,
  '2024-06-01', 101,
  '2024-07-01', 122,
  '2024-08-01', 136,
  '2024-09-01', 156,
  '2024-10-01', 175,
  '2024-11-01', 204,
  '2024-12-01', 224,
  '2025-01-01', 243, # new year
]

# 2025 monthly publication totals (first of each month)
pub_2025 = [
  '2024-09-01', 0,
  '2024-10-01', 4,
  '2024-11-01', 8,
  '2024-12-01', 13,
  '2025-01-01', 21, # new year
]

"""
Function that plots one year's worth of data
"""
def plot_one_year(year, data, color):
    x = []
    y = []
    start_date = datetime.datetime.strptime(str(year) + '-01-01', '%Y-%m-%d')
    for i in range(0, len(data), 2):
        end_date = datetime.datetime.strptime(data[i], '%Y-%m-%d')
        # Compute the number of months between two dates:
        num_months = (end_date.year - start_date.year) * 12 + (end_date.month - start_date.month)
        x.append(num_months)
        y.append(data[i+1])

    # Plot data
    ax1.plot(x, y, marker='o', color=color, linewidth=2, label=str(year))

# Parse command line args
parser = argparse.ArgumentParser()
parser.add_argument("--png", action='store_true', default=False)
args = parser.parse_args()

# Draw vertical lines at beg/end of year. Note: the data is from
# the first of every month, so draw the end of year line on
# Jan 1st of the following year.
ax1.axvline(x=0, linewidth=1, linestyle='--', color='lightgrey')
ax1.axvline(x=12, linewidth=1, linestyle='--', color='lightgrey')

# Plot the monthly data for each year
plot_one_year(2017, pub_2017, muted_dark_blue)
plot_one_year(2018, pub_2018, muted_light_blue)
plot_one_year(2019, pub_2019, muted_green)
plot_one_year(2020, pub_2020, muted_red)
plot_one_year(2021, pub_2021, muted_pink)
plot_one_year(2022, pub_2022, 'k')
plot_one_year(2023, pub_2023, muted_grey)
plot_one_year(2024, pub_2024, muted_yellow)

# Label beginning and end of year.
ax1.set_xticks([0, 12])
ax1.set_xticklabels(['Year Start', 'Year End'])

# Axis labels
ax1.set_ylabel('N. Publications')
ax1.legend(loc='upper left', prop=fontP)

# Write to PDF unless script was run with --png option.
# On Ubuntu, you may need to run the following
# command to get this working:
# sudo apt-get install dvipng
# To subsequently update the website,
# cp *.png ~/projects/libMesh.github.io/images/
# and then push the changes.
if args.png:
    plt.savefig('libmesh_citations_monthly.png', format='png', dpi=200)
else:
    plt.savefig('libmesh_citations_monthly.pdf', format='pdf')
