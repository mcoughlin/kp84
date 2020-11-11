## Code to re-run Sextractor on the LCO images 
## Change the folder that contains the unpacked LCO data

import optparse
import warnings
import glob
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import os

from astropy.table import Table

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

__version__ = 1.0

def parse_commandline():
    """@Parse the options given on the command-line.
    """
    parser = optparse.OptionParser(usage=__doc__,version=__version__)

    parser.add_option("-l","--transients", help="Transient list.", default='../data/Jacquinery/transients.dat')

    parser.add_option("-v", "--verbose", action="store_true", default=False,
                      help="Run verbosely. (Default: False)")

    opts, args = parser.parse_args()

    # show parameters
    if opts.verbose:
        print >> sys.stderr, ""
        print >> sys.stderr, "running gwemopt_run..."
        print >> sys.stderr, "version: %s"%__version__
        print >> sys.stderr, ""
        print >> sys.stderr, "***************** PARAMETERS ********************"
        for o in opts.__dict__.items():
          print >> sys.stderr, o[0]+":"
          print >> sys.stderr, o[1]
        print >> sys.stderr, ""

    return opts

# =============================================================================
#
#                                    MAIN
#
# =============================================================================

warnings.filterwarnings("ignore")

# Parse command line
opts = parse_commandline()

names = ('name', 'ra', 'dec')
transients = Table.read(opts.transients, format='ascii',
                        names=names, data_start=0)

for transient in transients:
    name = transient["name"]

    system_command = "python sextractor_mag.py -t %s -d ../data/Jacquinery/%s/ --doAstrometryNet" % (name, name)
    os.system(system_command)

