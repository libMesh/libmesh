#!/usr/bin/env python
from github_traffic_plotviews import *
from github_traffic_plotclones import *

# Actually make the plots.
PlotViews().plot_data()
PlotClones().plot_data()

# Local Variables:
# python-indent: 2
# End:
