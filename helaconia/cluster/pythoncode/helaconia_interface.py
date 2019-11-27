################################################################################
# 
# Running interface for HELAC-Onia on a cluster or a multicore
#
################################################################################
"""A user friendly command line interface to access HELAC-Onia features.
   Uses the cmd package for command interpretation and tab completion.
"""
from __future__ import division

import atexit
import cmath
import cmd
import glob
import logging
import logging.config
import coloring_logging
import math
import optparse
import os
import pydoc
import random
import re
import shutil
import signal
import stat
import subprocess
import sys
import time
import traceback
import tarfile
import copy
import datetime
import itertools

try:
    import readline
    GNU_SPLITTING = ('GNU' in readline.__doc__)
except:
    GNU_SPLITTING = True

root_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
root_path = os.path.split(root_path)[0]
sys.path.insert(0, os.path.join(root_path,'bin'))

# usefull shortcut
pjoin = os.path.join
# Special logger for the Cmd Interface
logger = logging.getLogger('helaconia.stdout') # -> stdout
logger_stderr = logging.getLogger('helaconia.stderr') # ->stderr

import extended_cmd as cmd
import misc as misc
import files as files
import cluster as cluster
import helaconia_run_interface as run_interface
from __init__ import InvalidCmd, HELACOniaError, HODIR

#===============================================================================
# HELACOniaCmd
#===============================================================================
class HELACOniaCmd(cmd.Cmd):

    options = {"pt_list":[],
               "subprocesses":[], # if empty, uses default one
               "decays":[], # if empty, uses default one
               'combine_factors':[] # if empty, uses [1.,1.,...]
        }

    addon_process = {}

    _display_opts = ["option"]

    def __init__(self, *arg, **opt):
        """ Special init tasks for the Loop Interface """
        logging.config.fileConfig(os.path.join(root_path, 'cluster', 'pythoncode', \
                                                   '.ho_logging.conf'))
        logging.root.setLevel(eval('logging.INFO'))
        info = misc.get_pkg_info()
        info_line = ""
        if info.has_key('version') and  info.has_key('date'):
            len_version = len(info['version'])
            len_date = len(info['date'])
            if len_version + len_date < 30:
                info_line = "#=         VERSION %s %s %s         =\n" % \
                            (info['version'],
                            (30 - len_version - len_date) * ' ',
                            info['date'])
        if info_line:
            info_line = info_line[1:]
        logger.info(\
        "\n============================================================\n" + \
        "=                                                          =\n" + \
        "=                     H E L A C - O n i a                  =\n" + \
        "=                                                          =\n" + \
        info_line + \
        "=                                                          =\n" + \
        "=                                                          =\n" + \
        "=       Author:        Hua-Sheng Shao                      =\n" + \
        "=   Affliation:         CERN, PH-TH                        =\n" + \
        "=      Contact:      erdissshaw@gmail.com                  =\n" + \
        "=    Reference:   Comput.Phys.Commun. 184 (2013) 2562      =\n" + \
        "=        arXiv:           1212.5293                        =\n" + \
        "=    Reference:   Comput.Phys.Commun. 198 (2016) 238       =\n" + \
        "=        arXiv:           1507.03435                       =\n" + \
        "=    Webpage: http://helac-phegas.web.cern.ch/helac-phegas =\n" + \
        "=                                                          =\n" + \
        "============================================================")

        self.set_addon_process() # read addon information

        self.set_options() # read options from default.inp and user.inp

        cmd.Cmd.__init__(self, *arg, **opt)

    def preloop(self):
        """Initializing before starting the main loop"""
        self.prompt = 'HO>'
        # preloop mother
        cmd.Cmd.preloop(self)

    @staticmethod
    def split_arg(line):
        """Split a line of arguments"""
        
        split = line.split()
        out=[]
        tmp=''
        for data in split:
            if data[-1] == '\\':
                tmp += data[:-1]+' '
            elif tmp:
                tmp += data
                tmp = os.path.expanduser(os.path.expandvars(tmp))
                out.append(tmp)
            else:
                out.append(data)
        return out

    def check_set(self, args, log=True):
        """ check the validity of the line"""

        if len(args) > 1 and '=' == args[1]:
            if not args[0] in self.options:
                raise HELACOniaError("option %s is not a valid one !"%args[0])
            args.pop(1)
            return
        raise HELACOniaError("set is wrong !")

    def check_generate(self,args,log=True):
        """check the validity of the line"""

        if 'addon' in args:
            if len(args)<2:
                raise HELACOniaError("Please specify the addon process number !\nSee addon/addon_process.dat.")
            if not int(args[1]) in self.addon_process.keys():
                raise HELACOniaError("The addon number=%s is wrong !\nSee addon/addon_process.dat."%args[1])
            return

        if len(args) > 2 and ">" == args[2]:
            args.pop(2)
            return
        if len(args) > 1 and ">" == args[1]:
            raise HELACOniaError("HELAC-Onia is unable to handle decay process !")

        raise HELACOniaError("generate is wrong !")

    def do_generate(self,line):
        """generate a process via shell"""
        args = self.split_arg(line)

        self.check_generate(args)

        if 'addon' in args:
            addon_subproc = self.addon_process[int(args[1])]
            self.options["subprocesses"].append(addon_subproc)
            return
                                            

        replace_mp_ids = copy.copy(ho_multparticles_id)

        replace_ids = copy.copy(ho_particles_id)
        replace_ids.update(ho_charmonia_id)
        replace_ids.update(ho_bottomonia_id)
        replace_ids.update(ho_Bc_id)

        #subprocesses = []
        subprocs = []
        for part in args:
            if (not part.lower() in replace_ids) and (not part.lower() in replace_mp_ids):
                raise HELACOniaError("%s is unknown !"%part)
            if part.lower() in replace_ids:
                multiparts=[replace_ids[part.lower()]]
            elif part.lower() in replace_mp_ids:
                multiparts=replace_mp_ids[part.lower()]
            if not subprocs:
                subprocs=[[part] for part in multiparts]
            else:
                newsubprocs = []
                for subproc,newpart in itertools.product(subprocs,multiparts):
                    subsubproc=subproc+[newpart]
                    newsubprocs+=[subsubproc]
                subprocs = copy.copy(newsubprocs)
                
        #subprocesses.append(subproc)
        
        self.options["subprocesses"] += subprocs
        
        #self.options["subprocesses"].append(subproc)

        #self.options["subprocesses"]=subprocesses

    def check_display(self, args):
        """check the validity of line
        syntax: display XXXXX
        """
        synstr="""ERROR in display syntax.
Syntax:
display option pt_list"""
        if len(args) < 1:
            raise HELACOniaError(synstr)
        if args[0] not in self._display_opts:
            raise HELACOniaError("Display options are %s"%" ".join(self._display_opts))
        if args[0] == "option":
            if len(args)<2:
                raise HELACOniaError(synstr)
            if args[1] not in self.options:
                raise HELACOniaError("%s is not a proper option"%args[1])

    def do_display(self,line):
        """display the options"""
        args = self.split_arg(line)

        self.check_display(args)

        if args[0] == "option":
            logger.info("%s:%s"%(args[1],str(self.options[args[1]])))

        return

    def check_decay(self,args,log=True):
        """check the validity of the line"""
        synstr="""ERROR in decay syntax.
Syntax:
e.g. decay cc~(3S11) > m+ m- @ 0.06d0"""
        if (">" not in args) or ("@" not in args):
            raise HELACOniaError(synstr)

        at_pos = args.index("@")
        if len(args)<at_pos+2:
            raise HELACOniaError(synstr)

        decay_BR=args[at_pos+1]

        
        decay_proc = []
        if len(args) > 1 and ">" == args[1]:
            decay_proc.append(args[0])
            decay_proc.extend(args[2:at_pos])
        else:
            raise HELACOniaError(synstr)

        return [decay_proc,decay_BR]
            

    def do_decay(self,line):
        """set the decay card via shell"""
        args = self.split_arg(line)

        decay_proc, decay_BR = self.check_decay(args)

        replace_ids = copy.copy(ho_particles_id)
        replace_ids.update(ho_charmonia_id)
        replace_ids.update(ho_bottomonia_id)
        replace_ids.update(ho_Bc_id)

        subproc = []
        for part in decay_proc:
            if not part.lower() in replace_ids:
                raise HELACOniaError("%s is unknown !"%part)
            subproc.append(replace_ids[part.lower()])

        self.options["decays"].append([subproc,decay_BR])

    def do_quit(self,line):
        """Not in help: exit the mainloop() """
        logger.info('Thanks for using HELAC-Onia ! Bye !')
        return True

    do_EOF = do_quit
    do_exit = do_quit
    # Set an option
    def do_set(self, line, log=True):
        """Set an option, which will be default for coming generations/outputs
        """
        # Be carefull:
        args = self.split_arg(line)

        # Check the validity of the arguments
        self.check_set(args)

        if args[0] == "pt_list":
            if len(args)>1:
                self.options[args[0]] = [ None if pt.lower() == "none" else pt for pt in args[1:]] 
                logger.info('set pt_list = [%s]'%(','.join(args[1:])))
            else:
                self.options[args[0]] = []
                logger.info('set pt_list = []')
            return

        if args[0] == "combine_factors":
            if len(args)>1:
                self.options[args[0]] = [ float(factor) for factor in args[1:]]
                logger.info('set combine_factors = [%s]'%(','.join(args[1:])))
            else:
                self.options[args[0]] = []
                logger.info('set combine_factors = []')
            return
        
        self.options[args[0]] = args[1]
        logger.info('set %s = %s'%(args[0],args[1]))
        return

    def set_addon_process(self):
        """assign addon processes from ./addon/addon_process.dat."""
        
        HODIR = root_path
        config_path = pjoin(HODIR,'addon','addon_process.dat')
        config_file = open(config_path)
        
        # read the file and extract information
        logger.info('load addon information from %s '%config_file.name)
        for line in config_file:
            split = line.split()
            newsplit = copy.copy(split)
            for i,data in enumerate(split):
                if data[0] == "#":
                    del newsplit[i:]
                    break
            if len(newsplit) > 1:
                self.addon_process[int(newsplit[0])]=newsplit[1]

        config_file.close()
        return self.addon_process

    def set_options(self):
        """assign all default variables from file ./input/default.inp and ./input/user.inp."""
        
        HODIR = root_path
        config_path = pjoin(HODIR, 'input', 'default.inp')
        config_file = open(config_path)
        
        # read the file and extract information
        logger.info('load initial variables from %s '% config_file.name)
        for line in config_file:
            split = line.split()
            newsplit = copy.copy(split)
            for i,data in enumerate(split):
                if data[0] == "#":
                    del newsplit[i:]
                    break
            if len(newsplit) > 1:
                if not newsplit[0] == "ptdisQ" and not newsplit[0] == "Pt1":
                    self.options[newsplit[0]]=newsplit[1]

        config_file.close()

        config_path = pjoin(HODIR, 'input', 'user.inp')
        config_file = open(config_path)
        logger.info('update initial variables from %s '% config_file.name)
        for line in config_file:
            split = line.split()
            newsplit = copy.copy(split)
            for i,data in enumerate(split):
                if data[0] == "#":
                    del newsplit[i:]
                    break
            if len(newsplit) > 1:
                if newsplit[0] in self.options:
                    self.options[newsplit[0]]=newsplit[1]

        config_file.close()

        return self.options

    def do_launch(self, line):
        """launch command"""
        run_interface.HELACOniaRunCmd(None,self.options)

    def do_shower(self, line):
        """run the shower program on a/several given parton level LHEF"""
        # HO> shower PROC_HO_0
        argss= self.split_arg(line)
        if os.path.exists(pjoin(os.getcwd(),argss[0])):
            self.options['shower_path']=pjoin(os.getcwd(),argss[0])
            run_interface.HELACOniaRunCmd(None,self.options)

    def cmdloop(self, intro=None):
        cmd.Cmd.cmdloop(self, intro)
        return

ho_multparticles_id = {'cc~(3pj1)':[443101,443111,443121],'cc~(3pj8)':[443108,443118,443128],'bb~(3pj1)':[553101,553111,553121],'bb~(3pj8)':[553108,553118,553128],'cb~(3pj1)':[453101,453111,453121],'cb~(3pj8)':[453108,453118,453128],'bc~(3pj1)':[-453101,-453111,-453121],'bc~(3pj8)':[-453108,-453118,-453128],'j':[3,4,35,-3,-4,8,-8],'p':[3,4,35,-3,-4,8,-8],'p~':[3,4,35,-3,-4,8,-8]}

ho_particles_id = {'ve':1, 'e-':2, 'u':3,'d':4,'vm':5,'m-':6,'c':7,'s':8,'vt':9,'tt-':10,'t':11,'b':12,\
                  'a': 31, 'z':32, 'w+':33,'w-':34,'g':35,'h':41,'g0':42,'g+':43,'g-':44,'e+':-2,'m+':-6,\
                      'tt+':-10,'ve~':-1,'u~':-3,'d~':-4,'vm~':-5,'c~':-7,'s~':-8,'vt~':-9,'t~':-11,'b~':-12}

ho_charmonia_id = {'cc~(1s01)':441001,'cc~(1s08)':441008,'cc~(3s11)':443011,'cc~(3s18)':443018,\
                  'cc~(1p11)':441111,'cc~(1p18)':441118,'cc~(3p01)':443101,'cc~(3p08)':443108,\
                  'cc~(3p11)':443111,'cc~(3p18)':443118,'cc~(3p21)':443121,'cc~(3p28)':443128}

ho_bottomonia_id = {'bb~(1s01)':551001,'bb~(1s08)':551008,'bb~(3s11)':553011,'bb~(3s18)':553018,\
                        'bb~(1p11)':551111,'bb~(1p18)':551118,'bb~(3p01)':553101,'bb~(3p08)':553108,\
                        'bb~(3p11)':553111,'bb~(3p18)':553118,'bb~(3p21)':553121,'bb~(3p28)':553128}

ho_Bc_id = {'cb~(1s01)':451001,'cb~(1s08)':451008,'cb~(3s11)':453011,'cb~(3s18)':453018,\
                'cb~(1p11)':451111,'cb~(1p18)':451118,'cb~(3p01)':453101,'cb~(3p08)':453108,\
                'cb~(3p11)':453111,'cb~(3p18)':453118,'cb~(3p21)':453121,'cb~(3p28)':453128,\
                'bc~(1s01)':-451001,'bc~(1s08)':-451008,'bc~(3s11)':-453011,'bc~(3s18)':-453018,\
                'bc~(1p11)':-451111,'bc~(1p18)':-451118,'bc~(3p01)':-453101,'bc~(3p08)':-453108,\
                'bc~(3p11)':-453111,'bc~(3p18)':-453118,'bc~(3p21)':-453121,'bc~(3p28)':-453128}
