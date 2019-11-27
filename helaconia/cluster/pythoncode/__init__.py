################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
class HELACOniaError(Exception):
    """Exception raised if an exception is find 
    Those Types of error will stop nicely in the cmd interface"""

class InvalidCmd(HELACOniaError):
    """a class for the invalid syntax call"""

import os
import logging
import time

#Look for basic file position HODIR
HODIR = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                                                os.path.pardir,os.path.pardir))

if ' ' in HODIR:
   logging.critical('''\033[1;31mpath to HO: "%s" contains space. 
    This is likely to create code unstability. 
    Please consider changing the path location of the code\033[0m''' % HODIR)
   time.sleep(1)
ReadWrite = True

try:
    open(os.path.join(HODIR,'.test'),'w').write('test')
    os.remove(os.path.join(HODIR,'.test'))
except IOError:
    ReadWrite = False
