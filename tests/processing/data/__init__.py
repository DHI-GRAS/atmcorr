import os
import glob

here = os.path.abspath(os.path.dirname(__file__))

MTDFILES = {os.path.basename(p): p for p in glob.glob(os.path.join(here, '*.imd'))}
