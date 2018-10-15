#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math

# Import stuff for working with dates
from datetime import datetime
from matplotlib.dates import date2num, num2date
import calendar

# Github has a "traffic" page now, but it doesn't seem like you can
# put in an arbitrary date range?  When I looked at it, it showed the
# numbers of unique visitors and views for the last two weeks only...
# So as long as I check back at least once every two weeks I can keep
# a record of it?  Great...

# https://github.com/libMesh/libmesh/graphs/traffic

# A generic class which can be used to customize the axis labels,
# titles, etc. of the three column plot data.
class PlotData(object):
  def __init__(self):
    self.left_axis_label = ''
    """
    Axis which goes on the left-hand side
    """

    self.right_axis_label = ''
    """
    Axis which goes on the right-hand side
    """

    self.weekly_plot_filename = ''
    """
    Filename for weekly plot.
    """

    self.monthly_plot_filename = ''
    """
    Filename for monthly plot.
    """

    self.title_string1 = ''
    """
    Part of title string describing the total
    """

    self.title_string2 = ''
    """
    Part of title string describing the average
    """

    self.data_file = ''
    """
    Name of a CSV file with date, total, unique data.
    """

  # Function which plots the two datasets on a single plot with two y axes.
  def plot_data(self):
    # Read the data from the CSV file. "|S11" refers to an 11 character string.
    # When you specify the dtype argument, the data is read into a structured
    # array and the columns can be accessed with the 'f0', 'f1', etc. names.
    data = np.genfromtxt(self.data_file, delimiter=',', dtype=("|S11", int, int))
    date_strings = data['f0']
    data_column2 = data['f1']
    data_column3 = data['f2']

    # Convert date strings into numbers.
    date_nums = []
    for d in date_strings:
      date_nums.append(date2num(datetime.strptime(d, '%Y-%b-%d')))

    # Initialize an array with 1, 7, 14, ...
    N = len(date_strings)
    week_indexes = range(0, N, 7)

    # Get total views and average unique viewers for each week
    weekly_column2 = []
    weekly_column3 = []
    x_axis = []
    for i in range(0, len(week_indexes)-1):
      start = week_indexes[i]
      stop = week_indexes[i+1]
      weekly_column2.append(sum(data_column2[start:stop]));
      weekly_column3.append(np.mean(data_column3[start:stop]));
      x_axis.append(date_nums[week_indexes[i]])

    # Get a reference to the figure
    fig = plt.figure()

    # The colors used come from sns.color_palette("muted").as_hex() They
    # are the "same basic order of hues as the default matplotlib color
    # cycle but more attractive colors."
    muted_dark_blue = u'#4878cf'
    muted_green = u'#6acc65'

    # 111 is equivalent to Matlab's subplot(1,1,1) command
    ax1 = fig.add_subplot(111)
    ax1.plot(x_axis, weekly_column2, color=muted_dark_blue, marker='o', linestyle='-', markersize=3)
    ax1.set_ylabel(self.left_axis_label + ' (blue circles)')

    # Choose the number of labels to create, then use linspace to create
    # them and convert them to ints.
    n_labels = 4
    x_axis_ticks = np.linspace(0, len(x_axis)-1, n_labels).astype(int)

    # Set tick labels and positions
    ax1.set_xticks([x_axis[i] for i in x_axis_ticks])
    ax1.set_xticklabels(['Week of \n' + str(num2date(x_axis[i]).date()) for i in x_axis_ticks])

    # Plot the weekly column 3 data.
    ax2 = ax1.twinx()
    ax2.plot(x_axis, weekly_column3, color=muted_green, marker='s', linestyle='-', markersize=3)
    ax2.set_ylabel(self.right_axis_label + ' (green squares)')

    # Add title
    title_string = self.title_string1 \
                   + ' ' \
                   + str(sum(data_column2)) \
                   + ', ' \
                   + self.title_string2 \
                   + ' ' \
                   + '%.1f' % np.mean(data_column3)
    fig.suptitle(title_string)

    # Save as PDF
    plt.savefig(self.weekly_plot_filename)



    # Make monthly plot
    fig.clf()

    # Generate date numbers at monthly intervals starting from '2014-Feb-17'
    now = datetime.now()
    month_intervals = [735281] # date2num for '2014-Feb-17'
    for yr in xrange(2014, now.year+1):
      for mo in xrange(1, 13):
        # Skip Jan 2014
        if (yr==2014 and mo==1):
          continue

        # http://stackoverflow.com/questions/27814983/how-to-get-month-interval-using-datetime-in-python
        # Add the number of days in the current month to the previous date num to
        # get the next one.
        month_intervals.append(month_intervals[-1] + calendar.monthrange(yr, mo)[1])

    # Find the indexes of each date.
    month_indexes = []
    for date in month_intervals:
      # Not all data sets have all of these dates, so we just use the
      # ones we have.
      if date in date_nums:
        month_indexes.append(date_nums.index(date))

    # Get total views and average unique viewers for each month
    monthly_column2 = [];
    monthly_column3 = [];
    x_axis = []
    for i in range(0, len(month_indexes)-1):
      start = month_indexes[i]
      stop = month_indexes[i+1]
      monthly_column2.append(sum(data_column2[start:stop]));
      monthly_column3.append(np.mean(data_column3[start:stop]));
      x_axis.append(date_nums[month_indexes[i]])

    # 111 is equivalent to Matlab's subplot(1,1,1) command
    ax1 = fig.add_subplot(111)
    ax1.plot(x_axis, monthly_column2, color=muted_dark_blue, marker='o', linestyle='-')
    ax1.set_ylabel(self.left_axis_label + ' (blue circles)')

    # Place an x-axis tick mark every x_step months.  As we get more data,
    # we'll have to increase x_step.
    x_step = 8
    x_axis_ticks = range(0, len(x_axis), x_step)

    # Set tick labels and positions
    ax1.set_xticks([x_axis[i] for i in x_axis_ticks])
    ax1.set_xticklabels([num2date(x_axis[i]).strftime('%b\n%Y') for i in x_axis_ticks])

    # Plot the average weekly unique visitors
    ax2 = ax1.twinx()
    ax2.plot(x_axis, monthly_column3, color=muted_green, marker='s', linestyle='-')
    ax2.set_ylabel(self.right_axis_label + ' (green squares)')

    # Save as PDF
    plt.savefig(self.monthly_plot_filename)

# Local Variables:
# python-indent: 2
# End:
