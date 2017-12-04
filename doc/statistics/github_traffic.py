#!/usr/bin/env python
from github_traffic_base import *

# Derive from the PlotData to plot views.
class PlotViews(PlotData):
  def __init__(self):
    super(PlotData, self).__init__()
    self.left_axis_label = 'Weekly page views'
    self.right_axis_label = 'Avg. Daily Unique Visitors'
    self.weekly_plot_filename = 'weekly_github_traffic.pdf'
    self.monthly_plot_filename = 'monthly_github_traffic.pdf'
    self.title_string1 = 'Total Pageviews:'
    self.title_string2 = 'Avg. Daily Unique Visitors:'
    self.data_file = 'github_traffic_plotviews.csv'

# Derive from the PlotData to plot clones.
class PlotClones(PlotData):
  def __init__(self):
    super(PlotData, self).__init__()
    self.left_axis_label = 'Clones per Week'
    self.right_axis_label = 'Avg. Daily Unique Clones'
    self.weekly_plot_filename = 'weekly_github_clones.pdf'
    self.monthly_plot_filename = 'monthly_github_clones.pdf'
    self.title_string1 = 'Total Clones:'
    self.title_string2 = 'Avg. Daily Unique Cloners:'
    self.data_file = 'github_traffic_plotclones.csv'

# Actually make the plots.
PlotViews().plot_data()
PlotClones().plot_data()

# Local Variables:
# python-indent: 2
# End:
