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
from __init__ import InvalidCmd, HELACOniaError, HODIR

#===============================================================================
# HELACOniaRunCmd
#===============================================================================
class HELACOniaRunCmd(cmd.Cmd):

    debug_output = 'HO_debug'

    # The three options categories are treated on a different footage when a
    # set/save configuration occur. current value are kept in self.options
    options_configuration = {
                       'topdraw_path':'./td',
                       'exrootanalysis_path':'./ExRootAnalysis',
                       'root_path': None,
                       'lhapdf_path': None,
                       'pythia8_path':None,
                       'hwpp_path':None,
                       'hepmc_path':None,
                       'fastjet_path':None,
                       'fortran_compiler':None,
                       'cluster_type': 'pbs',
                       'cluster_status_update': (600, 30),
                       'cluster_nb_retry':1,
                       'cluster_retry_wait':300,
                       'run_mode':2,
                       'cluster_queue':None,
                       'cluster_time':None,
                       'cluster_node_core':None,
                       'cluster_memory':None,
                       'nb_core': None,
                       'cluster_temp_path':None
                       }


    options_ho = {'pt_list':[],# pt_list: [] means total cross section; 
                               # [None] means using the default one
                  'subprocesses':[], # subprocesses: [] means default subprocess
                  'decays':{}, # decays: [] means default decay card
                  'combine_factors':[],
                  'stdout_level':None}

    default_options = {}


    def __init__(self, ho_dir, options, *args, **opts):
        """Init history and line continuation"""
        logging.config.fileConfig(os.path.join(root_path, 'cluster', 'pythoncode', \
                                                   '.ho_logging.conf'))
        logging.root.setLevel(eval('logging.INFO'))
        cmd.Cmd.__init__(self, *args, **opts)
        # Define current HELAC-Onia directory
        if ho_dir is None:
            ho_dir = root_path

        self.ho_dir = ho_dir
        self.options = options

        # usefull shortcut
        self.status = pjoin(self.ho_dir, 'status')
        self.error =  pjoin(self.ho_dir, 'error')
        self.ptdir = {}
        self.job_core = {}

        # Check that the directory is not currently running
        pid = os.getpid()

        #self.to_store = []
        self.run_name = None
        self.gener_dict = {}
        #self.run_tag = None
        #self.banner = None
        # Load the configuration file
        self.set_configuration()
        self.configure_run_mode(self.options['run_mode'])
        # Load the default.inp and user.inp
        self.set_options()

        # Get number of initial states
        #nexternal = open(pjoin(self.ho_dir,'input','process.inp')).read()
        #nexternal = int(nexternal.split()[0])
        #self.ninitial = 2

        keyboardstop=False
        if 'shower_path' in self.options and self.options['shower_path']:
            cwd=self.options['shower_path']
        else:
            # create the directories at the current directory
            cwd=self.create_directories()


            # run the program
            try:
                self.run("HO",{},cwd=cwd)
            except KeyboardInterrupt:
                self.stop_on_keyboard_stop()
                logger.info('Stop by user')
                keyboardstop=True
                pass
            else:
                # collect the results
                logger.info('Collecting results')
                result_path = self.collect_results(cwd=cwd)
                logger.info('Results are collected in %s'%result_path)

        if not keyboardstop:
            try:
                # run the shower
                self.run_lhe_shower("HO_Shower",{},cwd=cwd)
            except KeyboardInterrupt:
                keyboardstop=True
                pass

    def stop_on_keyboard_stop(self):
        """action to perform to close nicely on a keyboard interupt"""
        try:
            if hasattr(self, 'cluster'):
                logger.info('rm jobs on queue')
                self.cluster.remove()
        except:
            pass

    def collect_results(self,cwd=None):
        """collect the results for parallel computations."""
        p_dirs = [file for file in os.listdir(cwd)
                  if file.startswith('P') and \
                      os.path.isdir(pjoin(cwd,file))]
        result_dir = pjoin(cwd,"results")
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        ptQ=False
        if not self.options.has_key('pt_list') or len(self.options['pt_list']) == 0:
            # total cross section
            text = "# Total Cross Section"
            text += "\n# SubDir            NMAX    ACC                       SD                      "
        else:
            # pt distribution
            text = "# Pt Differential Distribution"
            text += "\n# SubDir            NMAX    Pt          ACC                       SD                      "
            ptQ = True

        for i,p_dir in enumerate(p_dirs):
            gener = self.gener_dict[p_dir]
            tmp_dir = p_dir.strip()
            tmp_len = len(tmp_dir)
            np_dir = " "*max(17,tmp_len)
            np_dir = tmp_dir+np_dir[tmp_len:]
            screenoutput=pjoin(cwd,p_dir,"output","screen_output.txt")
            tmpnmax, tmpacc, tmpsd =self.read_result(screenoutput,gener=gener)
            tmpnmax=tmpnmax.lower().strip()
            tmpacc=tmpacc.lower().strip()
            tmpsd =tmpsd.lower().strip()
            nmax = " "*max(7,len(tmpnmax))
            acc = " "*max(25,len(tmpacc))
            sd = " "*max(25,len(tmpsd))
            nmax =tmpnmax+nmax[len(tmpnmax):]
            acc =tmpacc+acc[len(tmpacc):]
            sd =tmpsd+sd[len(tmpsd):]
            stradd = [np_dir,nmax]
            if ptQ:
                if p_dir in self.ptdir:
                    tmppt=self.ptdir[p_dir].lower().replace("d","e").strip()
                else:
                    tmppt=self.options['pt_list'][i].lower().replace("d","e").strip()
                tmppt_len = len(tmppt)
                pt = " "*max(11,tmppt_len)
                pt =tmppt+pt[tmppt_len:]
                stradd += [pt]
            stradd += [ acc, sd ]
            text +="\n  "+" ".join(stradd)
        result_path = pjoin(result_dir,"results.out")
        result_file = open(result_path,"w")
        result_file.write(text)
        result_file.close()

        if not ptQ:
            self.collect_plots(cwd=cwd)

        return result_path

    def collect_plots(self,cwd=None,dirs=None,appstr=""):
        """collect the plots for parallel computations."""
        if dirs == []:return
        if not dirs:
            p_dirs = [file for file in os.listdir(cwd)
                      if file.startswith('P') and \
                          os.path.isdir(pjoin(cwd,file))]
        else:
            p_dirs = dirs
        result_dir = pjoin(cwd,"results")
        if not os.path.exists(result_dir):
            os.mkdir(result_dir)
        if not p_dirs:return
        combine_factors = self.options['combine_factors']
        td_plots = []
        gnu_plots = []
        root_plots = []
        for i,p_dir in enumerate(p_dirs):
            if not dirs:
                out_dir = pjoin(cwd,p_dir,"output")
            else:
                out_dir = p_dir
            td_plots.extend([pjoin(out_dir,file) \
                                 for file in os.listdir(out_dir)\
                                 if file.endswith('.top') and \
                                 os.path.isfile(pjoin(out_dir,file))])
            gnu_plots.extend([pjoin(out_dir,file) \
                                 for file in os.listdir(out_dir)\
                                 if file.endswith('.gnu') and \
                                 os.path.isfile(pjoin(out_dir,file))])
            root_plots.extend([pjoin(out_dir,file) \
                                 for file in os.listdir(out_dir)\
                                 if file.endswith('.C') and \
                                 os.path.isfile(pjoin(out_dir,file))])

        # combine topdrawer file
        if td_plots:
            self.combine_td_gnu(td_plots,result_dir,type="td",appstr=appstr)
        if gnu_plots:
            self.combine_td_gnu(gnu_plots,result_dir,type="gnu",appstr=appstr)
        if root_plots:
            self.combine_root(root_plots,result_dir,appstr=appstr)

    def combine_root(self,root_plots,result_dir,appstr=""):
        """combine root files"""
        rt_text = ""
        rt_combine_factors = self.options['combine_factors']
        if len(rt_combine_factors) < len(root_plots):
            rt_combine_factors.extend([1. for ij in \
                                           range(len(root_plots)-len(rt_combine_factors))])
        elif len(rt_combine_factors) > len(root_plots):
            del rt_combine_factors[len(root_plots):]

        store_plots_str = {}
        store_plots = {}
        store_index = {}

        dim_type=2
        dim_type_list = []
        for i_rt,rt_plot in enumerate(root_plots):
            rt_file = open(rt_plot)
            store_res = 0
            yfactor = rt_combine_factors[i_rt]
            i = 0
            if i_rt == 0:
                updata = []
                for line in rt_file:
                    if self.plot_start(line,type="root"):
                        store_res = 1
                        i = i + 1
                    elif self.plot_end(line,type="root"):
                        store_res = 0
                        # following are stable for empty plot 
                        if 'hist_3d -> SetOption("colz");' in line:
                            dim_type=3
                        else:
                            dim_type=2
                        if len(dim_type_list)<i:
                            dim_type_list.append(dim_type)

                    if store_res == 0:
                        if '\n' in line:
                            rt_text +='%s'%line
                        else:
                            rt_text +='%s\n'%line
                    elif self.plot_start(line,type="root"):
                        rep_str = 'plot_%i'%i
                        if '\n' in line:
                            rt_text +='%s'%line
                        else:
                            rt_text +='%s\n'%line
                        rt_text +='%('+rep_str+')s'
                        if not rep_str in store_plots:
                            store_plots[rep_str]=[]
                        store_index[rep_str]=i
                    else:
                        rep_str = 'plot_%i'%i
                        #if not rep_str in store_plots:
                        #    store_plots[rep_str]=[]
                        if 'xybin' in line:
                            dim_type=3
                        else:
                            dim_type=2
                        if len(dim_type_list)<i:
                            dim_type_list.append(dim_type)
                        data = line.split()
                        data.pop()
                        if 'FindBin' in line:
                            if dim_type==2:
                                updata.append(float(data[-1]))
                            else:
                                updata.append([float(data[-3]),float(data[-1])])
                        if 'SetBinContent' in line:
                            updata.append(float(data[-1])*yfactor)
                        if 'SetBinError' in line:
                            updata.append(float(data[-1])*yfactor)
                            store_plots[rep_str].append(updata)
                            #store_index[rep_str]=i
                            updata=[]
            else:
                for line in rt_file:
                    if self.plot_start(line,type="root"):
                        store_res = 1
                        i = i + 1
                    elif self.plot_end(line,type="root"):
                        store_res = 0

                    if store_res == 1 and not self.plot_start(line,type="root"):
                        if 'xybin' in line:
                            dim_type=3
                        else:
                            dim_type=2
                        rep_str = 'plot_%i'%i
                        data = line.split()
                        data.pop()
                        if 'FindBin' in line:
                            if dim_type==2:
                                j=self.find_pos(float(data[-1]),store_plots[rep_str])
                            else:
                                j=self.find_pos_3d([float(data[-3]),float(data[-1])],store_plots[rep_str])
                        if 'SetBinContent' in line:
                            store_plots[rep_str][j][1]+=float(data[-1])*yfactor
                        if 'SetBinError' in line:
                            store_plots[rep_str][j][2]+=float(data[-1])*yfactor

            rt_file.close()
        for key,value in store_plots.items():
            if dim_type_list[store_index[key]-1]==2:
                store_plots_str[key]=self.root_print(value,i=store_index[key])
            else:
                store_plots_str[key]=self.root_print_3d(value,i=store_index[key])

        rt_text = rt_text%store_plots_str
        res_filename = "results"+appstr+".C"
        result_path = pjoin(result_dir,res_filename)
        result_file = open(result_path,"w")
        result_file.write(rt_text)
        result_file.close()                

    def find_pos(self,x0,value_list):
        j=0
        for x,y,dy in value_list:
            if abs(x0-x)/max(abs(x0),abs(x),1e-99) < 1e-5:
                return j
            if x0 < x: break
            j=j+1
        value_list.insert(j,[x0,0.,0.])

        return j

    def find_pos_3d(self,xy0,value_list):
        j=0
        x0=xy0[0]
        y0=xy0[1]
        for xy,z,dz in value_list:
            x=xy[0]
            y=xy[1]
            if abs(x0-x)/max(abs(x0),abs(x),1e-99) < 1e-5\
                    and abs(y0-y)/max(abs(y0),abs(y),1e-99)<1e-5:
                return j
            if abs(x0-x)/max(abs(x0),abs(x),1e-99) < 1e-5 and y0 < y: break
            if x0 < x : break
            j=j+1
        value_list.insert(j,[[x0,y0],0.,0.])

        return j

    def combine_td_gnu(self,td_plots,result_dir,type="td",appstr=""):
        """combine topdrawer/gnuplot files"""
        td_text = ""
        td_combine_factors = self.options['combine_factors']
        if len(td_combine_factors) < len(td_plots):
            td_combine_factors.extend([1. for ij in \
                                           range(len(td_plots)-len(td_combine_factors))])
        elif len(td_combine_factors) > len(td_plots):
            del td_combine_factors[len(td_plots):]

        store_plots_str = {}
        store_plots = {}

        dim_type=2
        dim_type_list=[]
        for i_td,td_plot in enumerate(td_plots):
            td_file = open(td_plot)
            store_res = 0
            yfactor = td_combine_factors[i_td]
            i = 0
            if i_td == 0:
                for line in td_file:
                    if self.plot_start(line,type=type):
                        store_res = 1
                        i = i + 1
                        dim_type_list.append(self.plot_dimtype(line,type=type))
                    elif self.plot_end(line,type=type):
                        store_res = 0

                        
                    if store_res == 0:
                        if '\n' in line:
                            td_text +='%s'%line
                        else:
                            td_text +='%s\n'%line
                    elif self.plot_start(line,type=type):
                        rep_str = 'plot_%i'%i
                        if '\n' in line:
                            td_text +='%s'%line
                        else:
                            td_text +='%s\n'%line
                        td_text +='%('+rep_str+')s'
                    else:
                        rep_str = 'plot_%i'%i
                        if not rep_str in store_plots:
                            store_plots[rep_str]=[]
                        data = line.split()
                        data = [float(val) for val in data]
                        if len(data)>2:
                            if dim_type_list[i-1]==2:
                                data[1] = data[1]*yfactor
                                data[2] = data[2]*yfactor
                                store_plots[rep_str].append(data)
                            elif type=="gnu":
                                data[2] = data[2]*yfactor
                                store_plots[rep_str].append(data)
                        elif dim_type_list[i-1] == 3 and type == "gnu":
                            store_plots[rep_str].append(None)

            else:
                j = 0
                for line in td_file:
                    if self.plot_start(line,type=type):
                        store_res = 1
                        i = i + 1
                    elif self.plot_end(line,type=type):
                        store_res = 0
                        j = 0

                    if store_res == 1 and not self.plot_start(line,type=type):
                        rep_str = 'plot_%i'%i
                        data = line.split()
                        data = [float(val) for val in data]
                        if len(data)>2:
                            if dim_type_list[i-1]==2:
                                data[1] = data[1]*yfactor
                                data[2] = data[2]*yfactor
                                store_plots[rep_str][j][1]+=data[1]
                                store_plots[rep_str][j][2]+=data[2]
                                j = j + 1
                            elif type=="gnu":
                                data[2] = data[2]*yfactor
                                store_plots[rep_str][j][2]+=data[2]
                                j = j + 1
                        elif dim_type_list[i-1]==3 and type=="gnu":
                            j = j + 1
            td_file.close()

        for key,value in store_plots.items():
            store_plots_str[key]=self.plot_print(value,type=type)

        if type=="gnu":
            td_text=td_text.replace('set format y "10^{%T}"','set format y "10^{%%T}"')
        td_text = td_text%store_plots_str
        res_filename = "results"+appstr
        if type == "td":
            res_filename += ".top"
        elif type == "gnu":
            res_filename += ".gnu"
        else:
            res_filename +=".unknown"
        result_path = pjoin(result_dir,res_filename)
        result_file = open(result_path,"w")
        result_file.write(td_text)
        result_file.close()

        return

    def plot_dimtype(self,line,type="td"):
        if type == "td":
            # only 2-dim for topdrawer                                                  
            return 2
        elif type == "gnu":
            if 'splot "-" with pm3d' in line:
                return 3
            else:
                return 2

    def plot_start(self,line,type="td"):
        if type == "td":
            return self.td_start(line)
        elif type == "gnu":
            return self.gnu_start(line)
        elif type == "root":
            return self.root_start(line)
        else:
            return False

    def plot_end(self,line,type="td"):
        if type == "td":
            return self.td_end(line)
        elif type == "gnu":
            return self.gnu_end(line)
        elif type == "root":
            return self.root_end(line)
        else:
            return False

    def root_start(self,line):
        if 'SetStats(false)' in line:
            return True
        else:
            return False

    def root_end(self,line):
        if 'histos -> Add(hist);' in line or 'hist_3d -> SetOption("colz");' in line:
            return True
        else:
            return False

    def gnu_start(self,line):
        if 'plot "-" with histeps' in line or 'splot "-" with pm3d' in line:
            return True
        else:
            return False

    def gnu_end(self,line):
        newline=line.split()
        if len(newline) == 1 and newline[0] == 'e':
            return True
        else:
            return False

    def td_start(self,line):
        if 'SET ORDER X Y DY' in line:
            return True
        else:
            return False

    def td_end(self,line):
        if 'HISTO' in line:
            return True
        else:
            return False

    def plot_print(self,value_list,type="td"):
        """output topdrawer/gnuplot file"""
        if type == "td":
            return self.td_print(value_list)
        elif type == "gnu":
            return self.gnu_print(value_list)

        return ""

    def td_print(self,value_list):
        """output topdrawer file"""
        text = ""
        for data in value_list:
            text += "       %5.4e       %5.4e       %5.4e\n"%tuple(data)

        return text

    def gnu_print(self,value_list):
        """output gnuplot file"""
        text=""
        for data in value_list:
            if data==None:
                text +="     \n"
            else:
                text +="     %7.6e       %7.6e    %7.6e\n"%tuple(data)

        return text

    def root_print(self,value_list,i=1):
        """output root file"""
        text=""
        for data in value_list:
            text +='   int xbin = id%i->FindBin(   %18.17e      );\n'%(i,data[0])
            text +='  id%i -> SetBinContent( xbin,   %18.17e );\n'%(i,data[1])
            text +='  id%i -> SetBinError( xbin,    %18.17e      );\n'%(i,data[2])

        return text

    def root_print_3d(self,value_list,i=1):
        """output root file"""
        text=""
        for data in value_list:
            text +='   int xybin = id3d%i->FindBin(   %17.16e      ,   %17.16e      );\n'%(i,data[0][0],data[0][1])
            text +='  id3d%i -> SetBinContent( xybin,   %17.16e );\n'%(i,data[1])
            text +='  id3d%i -> SetBinError( xybin,    %17.16e      );\n'%(i,data[2])

        return text

    def split_jobs(self):
        """split jobs for running on the cluster"""
        pass

    def read_result(self,file,gener=1):
        """read the result from file.
        gener:0 phegas; 1, vegas."""
        if not isinstance(file,str) or not os.path.isfile(file):
            raise HELACOniaError("%s is not a valid file"%str(file))
        result_file = open(file)
        res_list = []
        preline=""
        for line in result_file:
            newline = line.strip()
            if gener == 1:
                # VEGAS output
                if len(newline) > 53:
                    if newline[:25] == "="*20+"NCALL":
                        ncal=newline[25:].strip("=")
                        if "UNWEIGHTING" in preline:
                            ncal=ncal+"UNW"
                        ncal_save=ncal
                if "ITERATION" in newline:
                    if newline[:9] == "ITERATION":
                        it=newline[9:].replace(":","").strip()
                        ncal=ncal_save +"IT"+it
                if "+\-" in newline:
                    acc,sd = newline.split("+\-")
                    acc = acc.strip()
                    sd = sd.strip()
                    res_list.append((ncal,acc,sd))
            elif gener == 0:
                # PHEGAS output
                if "sigma=" in newline:
                    acc,sdper,nunwei,nused,ntot = newline[7:].split()
                    res_list.append((ntot,acc,str(float(sdper.replace('D','E'))*float(acc.replace('D','E')))))
            else:
                raise HELACOniaError("unknown type of gener = %d"%gener)
            preline = line
        result_file.close()
        res_list.reverse()
        res_nmax = "nan"
        res_acc = "nan"
        res_sd = "nan"
        if len(res_list) == 0:
            logger.warning("no result in %s"%file)
        for ncal,acc,sd in res_list:
            res_nmax = ncal
            res_acc = acc
            res_sd = sd
            if not "nan" in acc.lower() and not "nan" in sd.lower():
                break
        else:
            logger.warning("no valid result in %s"%file)
        return res_nmax,res_acc,res_sd
                
                    
    ############################################################################
    def create_directories(self):
        """create directories for parallel computations."""
        # create the directory at the working directory
        cwd = os.getcwd()
        #create_dir = os.environ["PWD"]
        create_dir = cwd
        for i in range(0,1000000):
            tmp_dir = pjoin(create_dir,"PROC_HO_%d"%i)
            if not os.path.exists(tmp_dir) or i > 1000000:
                create_dir = tmp_dir
                break
        if not os.path.exists(create_dir):
            os.mkdir(create_dir)
        # cd to the root_path/cluster
        createsubdir = pjoin(root_path,"cluster")
        os.chdir(createsubdir)
        # create the subdirectories
        if not self.options.has_key('pt_list') or len(self.options['pt_list']) == 0:
            # total cross section
            #subdir = pjoin(create_dir,"P0_calc_0")
            #subprocess.call(["bash","create_subdir.sh",\
            #                         subdir])
            replacedict={"ptdisQ":"F"}
            for key,value in self.options.items():
                if key in self.default_options and not value == self.default_options[key]:
                    replacedict[key]=value
            #self.change_user_card(subdir,replacedict)
            if self.options.has_key('subprocesses')and self.options['subprocesses']:
                i = 0
                for subproc in self.options['subprocesses']:
                    # subproc = self.options['subprocesses'][0]
                    if isinstance(subproc,str):
                        subdir = pjoin(create_dir,"P0_addon_%s"%subproc)
                        addonexefile = self.link_addon(subproc)
                        subprocess.call(["bash","create_subdir.sh",\
                                             subdir,addonexefile])
                        # copy the special files in the addon 
                        if os.path.exists(pjoin(HODIR,'addon',subproc,'input')):
                            cwdcwd=os.getcwd()
                            os.chdir(pjoin(HODIR,'addon',subproc,'input'))
                            for file in glob.glob("*"):
                                shutil.copy(file,pjoin(subdir,'input',file))
                    else:
                        subdir = pjoin(create_dir,"P0_calc_%d"%i)
                        subprocess.call(["bash","create_subdir.sh",\
                                             subdir])
                    self.change_user_card(subdir,replacedict)
                    if not isinstance(subproc,str):
                        self.change_process_card(subdir,subproc)
                    self.change_decay_card(subdir,self.options["decays"])
                    i = i+1
            else:
                subdir = pjoin(create_dir,"P0_calc_0")
                subprocess.call(["bash","create_subdir.sh",\
                                     subdir])
                self.change_user_card(subdir,replacedict)
                self.change_decay_card(subdir,self.options["decays"])
        else:
            # pT distributions
            i = 0
            replacedict={"ptdisQ":"T","Pt1":"10d0"}
            for key,value in self.options.items():
                if key in self.default_options and not value ==self.default_options[key]:
                    replacedict[key]=value

            if self.options.has_key('subprocesses') and self.options['subprocesses']:
                if len(self.options['subprocesses'])>1:
                    logger.warning('Only calculate the first subprocess for Pt distribution')
                subproc = self.options['subprocesses'][0]
            else:
                subproc = None

            for pt in self.options['pt_list']:
                subdir = pjoin(create_dir,"P0_calc_pT_%d"%i)
                if isinstance(subproc,str):
                    addonexefile = self.link_addon(subproc)
                    subprocess.call(["bash","create_subdir.sh",\
                                         subdir,addonexefile])
                else:
                    subprocess.call(["bash","create_subdir.sh",\
                                         subdir])
                self.ptdir["P0_calc_pT_%d"%i] = pt
                replacedict["Pt1"]=pt
                self.change_user_card(subdir,replacedict)
                #if self.options.has_key('subprocesses') and self.options['subprocesses']:
                #    if len(self.options['subprocesses'])>1:
                #        logger.warning('Only calculate the first subprocess for Pt distribution')
                #    subproc = self.options['subprocesses'][0]
                if isinstance(subproc,str):
                    self.link_addon(subproc)
                elif not subproc:
                    self.change_process_card(subdir,subproc)
                self.change_decay_card(subdir,self.options["decays"])
                i = i+1

            
        # return to the working directory
        os.chdir(cwd)
        return create_dir

    def link_addon(self,subproc):
        """link the addon process subproc"""
        addonexefile="HO_%s"%subproc
        HODIR = root_path
        config_path = pjoin(HODIR,'bin',addonexefile)
        program = misc.which(config_path)
        if program:
            return addonexefile
        makefile="makefile_%s"%subproc
        config_path = pjoin(HODIR,'addon',subproc,makefile)
        if not os.path.isfile(config_path):
            raise HELACOniaError("Cannot find makefile in addon/%s"%subproc)
        shutil.copy(config_path,pjoin(HODIR,makefile))
        cwd = os.getcwd()
        os.chdir(HODIR)
        misc.compile(arg=['-f',makefile],cwd=HODIR,job_specs=False)
        #subprocess.call(['make','-f',makefile],cwd=HODIR)
        os.remove(pjoin(HODIR,makefile))
        os.chdir(cwd)
        config_path = pjoin(HODIR,'bin',addonexefile)
        program = misc.which(config_path)
        if program:
            return addonexefile
        else:
            raise HELACOniaError("Failed to compile for addon process %s"%subproc)

    def change_process_card(self,subdir,subproc):
        """change the information from process.inp card."""

        if not subproc:
            return

        if not os.path.exists(subdir):
            raise HELACOniaError('%s is not a valid path'%subdir)

        config_path=pjoin(subdir,"input","process.inp")
        if os.path.isfile(config_path):
            os.remove(config_path)

        text = "%d\n%s\n"%(len(subproc),' '.join([str(part_id) for part_id in subproc]))
        config_file=open(config_path,"w")
        config_file.write(text)
        config_file.close()

    def change_decay_card(self,subdir,decay_list):
        """change the information from decay_user.inp card."""
        if not decay_list:
            return
        
        if not os.path.exists(subdir):
            raise HELACOniaError('%s is not a valid path'%subdir)

        config_path=pjoin(subdir,"input","decay_user.inp")
        if os.path.isfile(config_path):
            os.remove(config_path)

        text="# DECAY CHANNELS\n# The syntax is\n# Decay Chian I\n# k BR\n# M D1 D2 ... Dk"
        decay_i = 0
        for decay_proc,decay_br in decay_list:
            decay_i +=1
            text+="\nDecay Chain %i"%decay_i
            text+="\n%i %s"%(len(decay_proc)-1,decay_br)
            text+="\n%s"%' '.join([str(part_id) for part_id in decay_proc])
        
        config_file=open(config_path,"w")
        config_file.write(text)
        config_file.close()

    def change_user_card(self,subdir,replacedict):
        """change the information from user.inp card."""

        #if 'lhapdf' in replacedict:
        #    if replacedict['lhapdf'] and replacedict['lhapdf'].lower()!='f':
        #        replacedict['lhapdf']='T'
        #    else:
        #        replacedict['lhapdf']='F'
        replacekeys = replacedict.keys()
        if not replacekeys:
            return

        if not os.path.exists(subdir):
            raise HELACOniaError('%s is not a valid path'%subdir)

        config_path=pjoin(subdir,"input","user.inp")

        if not os.path.isfile(config_path):
            raise HELACOniaError('File %s does not exist'%config_path)

        config_file = open(config_path)
        # read the file and extract information
        logger.info('load user setup from %s ' % config_file.name)
        text=""
        keylen = [len(key) for key in replacekeys]
        used_keys = replacekeys
        find_gener = False
        for line in config_file:
            newline=line
            for i,key in enumerate(replacekeys):
                if len(newline) < keylen[i]:
                    continue
                if newline[:keylen[i]] != key:
                    continue
                used_keys[i]=None
                newline=key+" "+replacedict[key]
                break
            if "gener" in newline and not find_gener:
                if newline[:5] == "gener":
                    gener=int(newline[5:].split("#")[0].strip())
                    if gener < 3:
                        self.gener_dict[os.path.basename(subdir)]=0
                    else:
                        self.gener_dict[os.path.basename(subdir)]=1
                    find_gener=True
            if not "\n" in newline:
                newline=newline+"\n"
            text+=newline
        used_keys=list(set(used_keys))
        if not find_gener:
            self.gener_dict[os.path.basename(subdir)]=0
        if None in used_keys:
            used_keys.remove(None)
        if not used_keys:
            for key in used_keys:
                newline=key+" "+replacedict[key]
                text+=newline+"\n"
        config_file.close()
        os.remove(config_path)
        config_file=open(config_path,"w")
        config_file.write(text)
        config_file.close()

    

    ############################################################################
    def get_pdf_input_filename(self,cwd):
        """return the name of the file which is used by the pdfset"""

        if hasattr(self, 'pdffiles') and self.pdffiles:
            return self.pdffiles
        else:
            self.pdffiles = []
            for pdfgroup in ['cteq','mrs']:
                if pdfgroup == 'cteq':
                    for pdfset in ["cteq5l.tbl","cteq5m.tbl","cteq6d.tbl","cteq6l.tbl",\
                                       "cteq6l1.tbl","cteq6m.tbl"]:
                        self.pdffiles.append(pjoin(cwd,pdfgroup,pdfset))
                if pdfgroup == "mrs":
                    for pdfset in ["mrsb.dat","mrse.dat","mrst2002nlo.dat"]:
                        self.pdffiles.append(pjoin(cwd,pdfgroup,pdfset))
            return self.pdffiles
            #else:
                # possible when using lhapdf
            #    pdffile = subprocess.Popen('%s --pdfsets-path' % self.options['lhapdf'],
            #            shell = True, stdout = subprocess.PIPE).stdout.read().strip()
            #    self.pdffiles = self.pdffiles.append(pdffile)
            #    return self.pdffiles
                
      
    ############################################################################
    def configure_run_mode(self, run_mode):
        """change the way to submit job 0: single core, 1: cluster, 2: multicore"""

        self.cluster_mode = int(run_mode)
        
        if not hasattr(self,"nb_core"):self.nb_core=None
        
        if int(run_mode) == 2:
            if not self.nb_core:
                import multiprocessing
                if not self.options["nb_core"]:
                    self.nb_core = multiprocessing.cpu_count()
                else:
                    self.nb_core = int(self.options['nb_core'])
            nb_core =self.nb_core
        elif int(run_mode) == 0:
            nb_core = 1


        if int(run_mode) in [0, 2]:
            self.options["nb_core"]=nb_core
            self.cluster = cluster.MultiCore(
                             **self.options)
        
        if int(self.cluster_mode) == 1:
            opt = self.options
            cluster_name = opt['cluster_type']
            self.cluster = cluster.from_name[cluster_name](**opt)
        
        self.options["cluster_partition"] = \
            self.cluster.partition(self.options["cluster_time"])


    def update_status(self, status, level, force=True,
                      error=False, starttime = None, update_results=False,
                      print_log=True):
        """ update the index status """

        if print_log:
            if isinstance(status, str):
                if '<br>' not  in status:
                    logger.info(status)
            elif starttime:
                running_time = misc.format_timer(time.time()-starttime)
                logger.info(' Idle: %s,  Running: %s,  Completed: %s [ %s ]' % \
                           (status[0], status[1], status[2], running_time))
            else:
                logger.info(' Idle: %s,  Running: %s,  Completed: %s' % status[:3])

        #if update_results:
        #    self.results.update(status, level, error=error)


    ############################################################################
    def set_configuration(self, config_path=None, initdir=None):
        """ assign all configuration variable from file
            ./input/ho_configuration.txt. assign to default if not define """

        if not hasattr(self, 'options') or not self.options:
            self.options = dict(self.options_configuration)
            self.options.update(self.options_ho)

        for key in self.options_configuration.keys():
            if key not in self.options:
                self.options[key]=self.options_configuration[key]

        for key in self.options_ho.keys():
            if key not in self.options:
                self.options[key]=self.options_ho[key]

        if not config_path:
            if self.options.has_key('ho_path') and self.options['ho_path']:
                HODIR = self.options['ho_path']
            else:
                self.options['ho_path'] = root_path
                HODIR = root_path
            config_file = pjoin(HODIR, 'input', 'ho_configuration.txt')
            self.set_configuration(config_path=config_file, initdir=HODIR)
            return
            
        config_file = open(config_path)

        # read the file and extract information
        logger.info('load configuration from %s ' % config_file.name)
        for line in config_file:
            if '#' in line:
                line = line.split('#',1)[0]
            line = line.replace('\n','').replace('\r\n','')
            try:
                name, value = line.split('=')
            except ValueError:
                pass
            else:
                name = name.strip()
                value = value.strip()
                if name.endswith('_path'):
                    path = value
                    if os.path.isdir(path):
                        self.options[name] = os.path.realpath(path)
                        continue
                    if not initdir:
                        continue
                    path = pjoin(initdir, value)
                    if os.path.isdir(path):
                        self.options[name] = os.path.realpath(path)
                        continue
                elif name.lower() == "lhapdf":
                    #self.options['lhapdf_path'] = value
                    if value.lower() == "none":
                        self.options['lhapdf_path'] = None
                    else:
                        if os.path.exists(value):
                            self.options['lhapdf_path'] = \
                                os.path.split(os.path.dirname(\
                                    os.path.realpath( value )))[0]
                        else:
                            self.options['lhapdf_path'] = None
                elif name.lower() == "fastjet":
                    if value.lower() == "none":
                        self.options['fastjet_path'] = None
                    else:
                        if os.path.exists(value):
                            self.options['fastjet_path'] = \
                                os.path.split(os.path.dirname(\
                                    os.path.realpath( value )))[0]
                        else:
                            self.options['fastjet_path'] = None
                else:
                    self.options[name] = value
                    if value.lower() == "none":
                        self.options[name] = None


        # Treat each expected input
        # delphes/pythia/... path
        for key in self.options:
            # Final cross check for the path
            if key.endswith('path'):
                path = self.options[key]
                if path is None:
                    continue
                if os.path.isdir(path):
                    self.options[key] = os.path.realpath(path)
                    continue
                path = pjoin(self.ho_dir, self.options[key])
                if os.path.isdir(path):
                    self.options[key] = os.path.realpath(path)
                    continue
                elif self.options.has_key('ho_path') and self.options['ho_path']:
                    path = pjoin(self.options['ho_path'], self.options[key])
                    if os.path.isdir(path):
                        self.options[key] = os.path.realpath(path)
                        continue
                self.options[key] = None
            elif key.startswith('cluster') and key != 'cluster_status_update':
                if key in ('cluster_nb_retry','cluster_wait_retry'):
                    self.options[key] = int(self.options[key])
                if hasattr(self,'cluster'):
                    del self.cluster
                pass
            #elif key not in ['stdout_level']:
                # Default: try to set parameter
                #try:
                #    self.do_set("%s %s --no_save" % (key, self.options[key]), log=False)
                #except self.InvalidCmd:
            #    logger.warning("Option %s from config file not understood" \
            #                       % key)

        # Configure the way to open a file:
        #misc.open_file.configure(self.options)
        #self.configure_run_mode(self.options['run_mode'])


        return self.options

    def set_options(self):
        """assign all default variables from file ./input/default.inp and ./input/user.inp."""

        HODIR = root_path
        config_path = pjoin(HODIR, 'input', 'default.inp')
        config_file = open(config_path)

        # read the file and extract information                                  
        logger.info('load default variables from %s '% config_file.name)
        for line in config_file:
            split = line.split()
            newsplit = copy.copy(split)
            for i,data in enumerate(split):
                if data[0] == "#":
                    del newsplit[i:]
                    break
            if len(newsplit) > 1:
                if not newsplit[0] == "ptdisQ" and not newsplit[0] == "Pt1":
                    self.default_options[newsplit[0]]=newsplit[1]

        config_file.close()

        config_path = pjoin(HODIR, 'input', 'user.inp')
        config_file = open(config_path)
        logger.info('update default variables from %s '% config_file.name)
        for line in config_file:
            split = line.split()
            newsplit = copy.copy(split)
            for i,data in enumerate(split):
                if data[0] == "#":
                    del newsplit[i:]
                    break
            if len(newsplit) > 1:
                if newsplit[0] in self.default_options:
                    self.default_options[newsplit[0]]=newsplit[1]

        config_file.close()

        return self.default_options

    def run_exe(self, exe, args, run_type, cwd=None):
        """this basic function launch locally/on cluster exe with args as argument.
        """
        # first test that exe exists:                                                                                                                 
        execpath = None
        if cwd and os.path.exists(pjoin(cwd, exe)):
            execpath = pjoin(cwd, exe)
        elif not cwd and os.path.exists(exe):
            execpath = exe
        else:
            raise HELACOniaError('Cannot find executable %s in %s' \
                % (exe, os.getcwd()))
        # check that the executable has exec permissions                                                                                              
        if int(self.cluster_mode) == 1 and not os.access(execpath, os.X_OK):
            subprocess.call(['chmod', '+x', exe], cwd=cwd)
        # finally run it
        if int(self.cluster_mode) == 0:
            #this is for the serial run
            if run_type == "Showering events...":
                fscreenout = open(pjoin(cwd,'shower_screen_output,txt'),'w')
            else:
                fscreenout = open(pjoin(cwd,"output/screen_output.txt"),'w')
            #fscreenerr = open(pjoin(cwd,"output/screen_error.txt"),"w")
            misc.call(['./'+exe], cwd=cwd,stdout=fscreenout,stderr=subprocess.STDOUT)
            fscreenout.close()
            #fscreenerr.close()
            self.ijob += 1
            if run_type == "Showering events...":
                self.update_status((max([self.njobs - self.ijob - 1, 0]),
                                    min([1, self.njobs - self.ijob]),
                                    self.ijob, run_type), level='shower')
            else:
                self.update_status((max([self.njobs - self.ijob - 1, 0]),
                                    min([1, self.njobs - self.ijob]),
                                    self.ijob, run_type), level='parton')
        #elif 'reweight' in exe:
                #Find the correct PDF input file                                                                                                      
        #        input_files, output_files = [], []
        #        input_files.append(self.get_pdf_input_filename(cwd))
        #        input_files.append(pjoin(os.path.dirname(exe), os.path.pardir, 'reweight_xsec_events'))
        #        input_files.append(args[0])
        #        output_files.append('%s.rwgt' % os.path.basename(args[0]))
        #        output_files.append('reweight_xsec_events.output')
        #        output_files.append('scale_pdf_dependence.dat')

        #        return self.cluster.submit2(exe, args, cwd=cwd,
        #                         input_files=input_files, output_files=output_files)

        #this is for the cluster/multicore run                                                                                                        
        #elif 'ajob' in exe:
            # check if args is a list of string                                                                                                       
        #    if type(args[0]) == str:
        #        input_files, output_files, args = self.getIO_ajob(exe,cwd, args)
                #submitting                                                                                                                           
        #        self.cluster.submit2(exe, args, cwd=cwd,
        #                     input_files=input_files, output_files=output_files)

                # keep track of folders and arguments for splitted evt gen                                                                            
        #        if len(args) == 4 and '_' in output_files[-1]:
        #            self.split_folders[pjoin(cwd,output_files[-1])] = [exe] + args

        else:
            if int(self.cluster_mode) == 1:
                fscreenout = None
            else:
                if run_type == "Showering events...":
                    fscreenout = open(pjoin(cwd,'shower_screen_output.txt'),'w')
                else:
                    fscreenout = open(pjoin(cwd,"output/screen_output.txt"),'w')
            if exe in self.job_core:
                cl_args = self.cluster.submit_args(self.job_core[exe])
            else:
                cl_args = []
            return self.cluster.submit(exe, cl_args+self.options["cluster_partition"], cwd=cwd,stdout=fscreenout,stderr=subprocess.STDOUT)


    def check_event_files(self):
        """check the integrity of the event files after splitting, and resubmit 
        those which are not nicely terminated"""
        to_resubmit = []
        for dir in self.split_folders.keys():
            last_line = ''
            try:
                last_line = subprocess.Popen('tail -n1 %s ' % \
                    pjoin(dir, 'events.lhe'), \
                shell = True, stdout = subprocess.PIPE).stdout.read().strip()
            except IOError:
                pass

            if last_line != "</LesHouchesEvents>":
                to_resubmit.append(dir)

        self.njobs = 0
        if to_resubmit:
            run_type = 'Resubmitting broken jobs'
            logger.info('Some event files are broken, corresponding jobs will be resubmitted.')
            logger.debug('Resubmitting\n' + '\n'.join(to_resubmit) + '\n')
            for dir in to_resubmit:
                files.rm([dir])
                job = self.split_folders[dir][0]
                args = self.split_folders[dir][1:]
                run_type = 'monitor'
                cwd = os.path.split(dir)[0]
                self.run_exe(job, args, run_type, cwd=cwd )
                self.njobs +=1

            self.wait_for_complete(run_type)

    def run(self, mode,options, cwd=None):
        "run HELACOniaRunCmd."
        logger.info('Starting run')
        if int(self.cluster_mode) == 1:
            cluster_name = self.options['cluster_type']
            self.cluster = cluster.from_name[cluster_name](**self.options)
        if int(self.cluster_mode) == 2:
            try:
                import multiprocessing
                if not self.nb_core:
                    try:
                        self.nb_core = int(self.options['nb_core'])
                    except TypeError:
                        self.nb_core = multiprocessing.cpu_count()
                logger.info('Using %d cores' % self.nb_core)
            except ImportError:
                self.nb_core = 1
                logger.warning('Impossible to detect the number of cores => Using One.\n'+
                        'Use set nb_core X in order to set this number and be able to'+
                                               'run in multicore.')
            self.cluster = cluster.MultiCore(**self.options)

        
        # find and keep track of all the jobs
        folder_names = {'HO': [ 'output*' ]}
        job_dict = {}
        p_dirs = [file for file in os.listdir(cwd)
                  if file.startswith('P') and \
                      os.path.isdir(pjoin(cwd,file))]
        if int(self.cluster_mode) != 1:
            for dir in p_dirs:
                job_dict[dir] = [file for file in os.listdir(pjoin(cwd,dir))\
                                     if file.startswith('Helac-Onia')]
        else:
            # running on the cluster
            # split the jobs first
            if self.options["cluster_node_core"] == None or not isinstance(self.options["cluster_node_core"],str):
                code_core = self.options["cluster_node_core"]
            elif isinstance(self.options["cluster_node_core"],str):
                if self.options["cluster_node_core"].lower() == "none":
                    code_core = None
                else:
                    code_core = int(self.options["cluster_node_core"])

            if not isinstance(code_core,int):
                code_core = 1
            chunks = lambda l,n: [l[x:x+n] for x in xrange(0,len(l),n)]
            p_dirs_list = chunks(p_dirs,code_core)
            bin_dir = pjoin(cwd,"bin")
            if not os.path.exists(bin_dir):
                os.mkdir(bin_dir)
            current_cwd = os.getcwd()
            for i, dirs in enumerate(p_dirs_list):
                #jobfile = open("ajob%d"%i,"w")
                text = "#!/bin/bash\n"
                for dir in dirs:
                    text += "cd %s\n"%pjoin(cwd,dir)
                    text += "./Helac-Onia 1>./%s 2>&1 &\n"\
                        %pjoin("output","screen_output.txt")
                text += "cd %s"%current_cwd
                offset = 0
                while offset<100000:
                    jobpath = pjoin(bin_dir,"ajob%d"%(i+offset))
                    if not os.path.isfile(jobpath):
                        break
                    else:
                        offset += 1
                jobpath = pjoin(bin_dir,"ajob%d"%(i+offset))
                jobfile = open(jobpath,"w")
                jobfile.write(text)
                jobfile.close()
                self.job_core["ajob%d"%(i+offset)] = len(dirs) 
            job_dict["bin"]= \
                [ file for file in os.listdir(bin_dir)\
                      if file.startswith('ajob') ]
                          

        #ho_status = ['Computing cross section','Generating events']
        ho_status = ['Computing cross section' ]
        devnull = os.open(os.devnull, os.O_RDWR)
        if mode in ['HO']:
            for i,status in enumerate(ho_status):
                # check if need to split jobs
                # at least one channel must have enough events
                #try:
                #    nevents_unweighted = open(pjoin(self.ho_dir,
                #                                    'nevents_unweighted')).read().split('\n')
                #except IOError:
                #    nevents_unweighted = []
                nevents_unweighted = []
                nevt_job = 0
                split = i == 1 and \
                    nevt_job >0 and \
                    any([int(l.split()[1])> int(nevt_job) \
                             for l in nevents_unweighted if l])
                #if split:
                #    misc.call([pjoin(self.ho_dir, 'bin', 'internal', 'split_jobs.py')] + \
                #                   [self.run_card['nevt_job']],
                #              stdout = devnull,
                #              cwd = self.ho_dir)
                #    assert os.path.exists(pjoin(self.ho_dir, 
                #            'nevents_unweighted_splitted'))

                self.update_status(status,level='parton')
                self.run_all(job_dict, [['2','B','%d' %i ]],status,split_jobs = split,cwd=cwd)

                # check that split jobs are all correctly terminated
                #if split:
                #    self.check_event_files()

        if int(self.cluster_mode) == 1:
            # if cluster run, wait 10 sec so that event files are  transferred back
            self.update_status(\
            "Waiting while files are transferred back from the cluster nodes",\
                level='parton')
            time.sleep(10)
        #if split:
        #    files.cp(pjoin(self.ho_dir,'nevents_unweighted_splitted'),\
        #                 pjoin(self.ho_dir,'nevents_unweighted'))

    def wait_for_complete(self, run_type,cwd = None):
        """this function waits for jobs on cluster to complete their run."""

        starttime = time.time()
        #logger.info('     Waiting for submitted jobs to complete')
        update_status = lambda i, r, f: self.update_status((i, r, f, run_type), 
                      starttime=starttime, level='parton', update_results=False)
        try:
            self.cluster.wait(cwd, update_status)
        except:
            self.cluster.remove()
            raise

    def run_all(self,job_dict,arg_list,run_type='monitor',split_jobs = False,cwd=None):
        """runs the jobs in job_dict (organized as folder: [job_list]), with arguments args"""
        self.njobs = sum(len(jobs) for jobs in job_dict.values())*len(arg_list)
        njob_split = 0
        self.ijob = 0
        if int(self.cluster_mode) == 0:
            self.update_status((self.njobs-1,1,0,run_type),level="parton")

        # this is to keep track, if splitting evt generation, of the various
        # folders/args in order to resubmit the jobs if some of them fail
        self.split_folders = {}
        for args in arg_list:
            for Pdir, jobs in job_dict.items():
                for job in jobs:
                    if not split_jobs:
                        self.run_exe(job,args,run_type,cwd=pjoin(cwd,Pdir))
                    else:
                        for n in self.find_jobs_to_split(Pdir,job,args[1]):
                            self.run_exe(job,args + [n],run_type,cwd=pjoin(cwd,Pdir))
                            njob_split += 1
        if int(self.cluster_mode) == 2:
            time.sleep(1) # security to allow all jobs to be launched
        if njob_split > 0:
            self.njobs = njob_split

        self.wait_for_complete(run_type,cwd = cwd)

    def get_option(self,option):
        """get the option from self.options or self.default_options"""
        if not option in self.options:
            if not option in self.default_options:
                raise HELACOniaError,"Unknow of option = %s"%option
            else:
                return self.default_options[option]
        else:
            return self.options[option]

    def create_shower_input(self,evt_file,cwd):
        """create the shower input script that will used by shower script."""
        shower = self.get_option('parton_shower')
        if isinstance(shower,str):shower=int(shower)
        shower_path=""
        if shower == 0:
            return None,None
        elif shower == 1:
            shower_tag = "PYTHIA8"
            MCmass_path=os.path.join(root_path,"shower","PYTHIA8","MCmasses_PYTHIA8.inc")
            if not os.path.exists(self.get_option('pythia8_path')):
                raise HELACOniaError, 'Cannot find pythia8 in %s'%self.get_option('pythia8_path')
            shower_dir=os.path.join(cwd,"PYTHIA8")
            shower_path = "PY8PATH=%s\n"%self.get_option('pythia8_path')
            #if not os.path.exists(self.get_option('hepmc_path')):
            #    raise HELACOniaError, 'Cannot find HepMC in %s'%self.get_option('hepmc_path')
        elif shower == 2:
            shower_tag = "PYTHIA6PT"
            shower_dir=os.path.join(cwd,"PYTHIA6")
            MCmass_path=os.path.join(root_path,"shower","PYTHIA6","MCmasses_PYTHIA6PT.inc")
        elif shower == 3:
            shower_tag = "PYTHIA6Q"
            shower_dir=os.path.join(cwd,"PYTHIA6")
            MCmass_path=os.path.join(root_path,"shower","PYTHIA6","MCmasses_PYTHIA6Q.inc")
        elif shower == 4:
            shower_tag = "HERWIG6"
            shower_dir=os.path.join(cwd,"HERWIG6")
            MCmass_path=os.path.join(root_path,"shower","HERWIG6","MCmasses_HERWIG6.inc")
        elif shower == 5:
            shower_tag = "HERWIGPP"
            shower_dir=os.path.join(cwd,"HERWIGPP")
            MCmass_path=os.path.join(root_path,"shower","HERWIGPP","MCmasses_HERWIGPP.inc")
        else:
            raise HELACOniaError, "Unknown the shower option: %d"%shower
        # loading the Monte Carlo masses
        if not os.path.exists(MCmass_path):
            raise HELACOniaError, "Cannot find the MC mass file %s"%MCmass_path
        logger.info("loading Monte Carlo Masses for %s"%shower_tag)
        mcmass_dict = {}
        mcmass_file = open(MCmass_path)
        for line in mcmass_file:
            pdg=int(line.split("=")[0].split("(")[1].split(")")[0])
            mass=float(line.split("=")[1].replace('d','e'))
            mcmass_dict[pdg] = mass

        logger.info("create shower input file for %s"%shower_tag)
        lhapdf = self.get_option('lhapdf')
        pdf = self.get_option('pdf')
        # the maximus number of events
        #nevents = self.get_option('nmc')
        
        content = 'ShowerMode=%s\n'%shower_tag
        content += 'ShowerDirectory=%s\n'%shower_dir
        if shower == 1 and self.get_option('hepmc_path') and \
                os.path.exists(self.get_option('hepmc_path')):
            content += 'HepMCDirectory=%s\n'%self.get_option('hepmc_path')
        content += 'HODirectory=%s\n'%HODIR
        content += 'EVTFILE=%s\n'%evt_file
        content += 'TMASS=%s\n'%self.get_option('tmass')
        content += 'TWIDTH=%s\n'%self.get_option('twidth')
        content += 'ZMASS=%s\n'%self.get_option('zmass')
        content += 'ZWIDTH=%s\n'%self.get_option('zwidth')
        content += 'WMASS=%s\n'%self.get_option('wmass')
        content += 'WWIDTH=%s\n'%self.get_option('wwidth')
        content += 'HGGMASS=%s\n'%self.get_option('higmass')
        content += 'HGGWIDTH=%s\n'%self.get_option('higwidth')
        content += 'beammon1=%s\n'%self.get_option('energy_beam1')
        content += 'beammon2=%s\n'%self.get_option('energy_beam2')
        if int(self.get_option('colpar')) == 1:
            # pp
            beam1=1
            beam2=1
        elif int(self.get_option('colpar')) == 2:
            # ppbar
            beam1=1
            beam2=-1
        elif int(self.get_option('colpar')) == 3:
            # ee
            beam1=0
            beam2=0
        else:
            raise HELACOniaError, "Wrong collision type : %s"%self.get_option('colpar')
        content += 'BEAM1=%d\n'%beam1
        content += 'BEAM2=%d\n'%beam2
        content += 'DMASS=%s\n'% mcmass_dict[1]
        content += 'UMASS=%s\n'% mcmass_dict[2]
        content += 'SMASS=%s\n'% mcmass_dict[3]
        content += 'CMASS=%s\n'% mcmass_dict[4]
        content += 'BMASS=%s\n'% mcmass_dict[5]
        content += 'EMASS=%s\n'% mcmass_dict[11]
        content += 'MUMASS=%s\n'% mcmass_dict[13]
        content += 'TAUMASS=%s\n'% mcmass_dict[15]
        content += 'GMASS=%s\n'% mcmass_dict[21]
        content += 'EVENT_NORM=average\n'
        if lhapdf.lower() in ['t','true','.true.']:
            content += 'PDFCODE=%s\n'%pdf
            content += 'LHAPDFPATH=%s\n'%self.get_option('lhapdf_path')
        else:
            content += 'PDFCODE=0\n'
        if self.get_option('ktrw').lower() in ['t','true','.true.']: 
            content += 'MLM_MERGING=1\n'
        else:
            content += 'MLM_MERGING=0\n'
        content += shower_path
        import shower_card
        wk_dir = os.path.split( cwd )[0]
        shower_path=pjoin(wk_dir,'input','shower_card_user.inp')
        if not os.path.isfile(shower_path):
            shower_path=pjoin(HODIR,'input','shower_card_user.inp')
        showercard=shower_card.ShowerCard(shower_path)
        content +=showercard.write_card(shower_tag)
        content +='\n'
        shower_input = pjoin(cwd,'HELACOnia_shower.inp')
        open(shower_input,'w').write(content)

        return shower_tag,shower_input

    def run_one_lhe_shower(self, evt_file,verbose=0):
        """runs shower program on the generated event LHEF, to produce showered events"""
        if not os.path.exists(evt_file):
            raise HELACOniaError,"Cannot find the event file %s"%evt_file
        wk_dir=os.path.split(os.path.dirname(os.path.realpath( evt_file )))[0]
        shower_dir = pjoin(wk_dir,'shower')
        shower,shower_input=self.create_shower_input(evt_file,shower_dir)
        if shower == None or shower_input == None:
            return None,None
        if not os.path.exists(shower_dir):
            raise HELACOniaError,'Cannot find the shower directory %s'%shower_dir

        extrapaths = []
        if shower == 'HERWIGPP':
            extrapaths.append(pjoin(self.get_option('hepmc_path'), 'lib'))
        if 'LD_LIBRARY_PATH' in os.environ.keys():
            ldlibrarypath = os.environ['LD_LIBRARY_PATH']
        else:
            ldlibrarypath = ''
        ldsplit=ldlibrarypath.split(':')
        extrapaths = [path for path in extrapaths if not path in ldsplit]
        ldlibrarypath += ':'+':'.join(extrapaths)
        os.putenv('LD_LIBRARY_PATH',ldlibrarypath)

        if verbose == 1:
            self.update_status('Compiling source codes for %s...' % shower, level='shower')
        # create a new shower running directory
        for i in range(0,1000000):
            tmp_dir = pjoin(shower_dir,"HO_%s_%d"%(shower,i))
            if not os.path.exists(tmp_dir) or i > 1000000:
                run_shower_dir = tmp_dir
                break
        if not os.path.exists(run_shower_dir):
            os.mkdir(run_shower_dir)
        shower_log = pjoin(run_shower_dir,'shower.log')
        #current_dir = os.getcwd()
        #os.chdir(shower_dir)
        output_log = open(shower_log,'w')
        subprocess.call(['bash','RUN_HO_SHOWER'],stdout=output_log,\
                  stderr=output_log,cwd=shower_dir)
        #os.chdir(current_dir)
        exe = 'HO_%s_SHOWER'% shower
        if not os.path.exists(pjoin(shower_dir,exe)) and \
                not os.path.exists(pjoin(shower_dir,'Pythia8.exe')):
            print open(shower_log).read()
            raise HELACOniaError,'Compilation failed, check %s for details'% shower_log
        logger.info('                     ... done')
        
        logger.info('Running in %s' % run_shower_dir)
        files.mv(shower_input,run_shower_dir)

        if shower != 'PYTHIA8':
            files.mv(pjoin(shower_dir, exe), run_shower_dir)
        else:
            exe='Pythia8.exe'
            files.mv(pjoin(shower_dir, 'Pythia8_lhe.cmnd'), run_shower_dir)
            files.mv(pjoin(shower_dir, 'Pythia8.exe'), run_shower_dir)
            files.ln(pjoin(self.get_option('pythia8_path'), 'examples', 'config.sh'),run_shower_dir)
            files.ln(pjoin(self.get_option('pythia8_path'),'xmldoc'),run_shower_dir)
        
        if shower == "HERWIGPP":
            try:
                files.ln(pjoin(self.get_option('hwpp_path'),'bin','Herwig++'),run_shower_dir)
            except Exception:
                raise HELACOniaError, 'The Herwig++ path set in the configuration file is not valid.'

        return run_shower_dir,exe

    def run_lhe_shower(self, mode,options, cwd=None):
        """run parton shower to generate showered events"""
        shower = self.get_option('parton_shower')
        if isinstance(shower,str):shower=int(shower)
        if shower == 0: return
        if shower == 1:
            shower_tag = "PYTHIA8"
        elif shower == 2:
            shower_tag = "PYTHIA6PT"
        elif shower == 3:
            shower_tag = "PYTHIA6Q"
        elif shower == 4:
            shower_tag = "HERWIG6"
        elif shower == 5:
            shower_tag = "HERWIGPP"
        unweight = self.get_option('unwgt')
        if unweight.lower() in ['f','false','.false.']: return
        # find all parton level LHEF
        logger.info("Preparing parton shower run")
        if int(self.cluster_mode) == 1:
            cluster_name = self.options['cluster_type']
            self.cluster = cluster.from_name[cluster_name](**self.options)
        if int(self.cluster_mode) == 2:
            try:
                import multiprocessing
                if not self.nb_core:
                    try:
                        self.nb_core = int(self.options['nb_core'])
                    except TypeError:
                        self.nb_core = multiprocessing.cpu_count()
                logger.info('Using %d cores' % self.nb_core)
            except ImportError:
                self.nb_core = 1
                logger.warning('Impossible to detect the number of cores => Using One.\n'+
                        'Use set nb_core X in order to set this number and be able to'+
                                               'run in multicore.')
            self.cluster = cluster.MultiCore(**self.options)

        # find and keep track of all the jobs
        #folder_names = {'HO': [ 'output*' ]}
        job_dict = {}
        run_shower_dirs = []
        run_shower_dict = {}
        p_dirs = [file for file in os.listdir(cwd)
                  if file.startswith('P') and \
                      os.path.isdir(pjoin(cwd,file))]
        
        for i,dir in enumerate(p_dirs):
            lhefile = [file for file in os.listdir(pjoin(cwd,dir,'output')) \
                           if file.endswith('.lhe')]
            if lhefile:
                lhefile=lhefile[0]
                evt_file=pjoin(cwd,dir,'output',lhefile)
                if i == 0:
                    verbose = 1
                else:
                    verbose = 0
                run_shower_dir,shower_exe = self.run_one_lhe_shower(evt_file,verbose=verbose)
                if run_shower_dir != None and shower_exe != None and os.path.exists(pjoin(run_shower_dir,shower_exe)):
                    run_shower_dict[run_shower_dir] = [ shower_exe ]
                    run_shower_dirs.append(run_shower_dir)
            else:
                if i == 0:
                    self.update_status('Compiling shower source codes ...', level='shower')

        if not run_shower_dict:
            return

        if not int(self.cluster_mode) == 1:
            job_dict = run_shower_dict
        else:
            # running on the cluster
            # split the jobs first
            if self.options["cluster_node_core"] == None or not isinstance(self.options["cluster_node_core"],str):
                code_core = self.options["cluster_node_core"]
            elif isinstance(self.options["cluster_node_core"],str):
                if self.options["cluster_node_core"].lower() == "none":
                    code_core = None
                else:
                    code_core = int(self.options["cluster_node_core"])

            if not isinstance(code_core,int):
                code_core = 1
            chunks = lambda l,n: [l[x:x+n] for x in xrange(0,len(l),n)]
            p_dirs_list = chunks(run_shower_dirs,code_core)
            bin_dir = pjoin(cwd,"bin")
            if not os.path.exists(bin_dir):
                os.mkdir(bin_dir)
            current_cwd = os.getcwd()
            for i, dirs in enumerate(p_dirs_list):
                #jobfile = open("ajob%d"%i,"w")
                text = "#!/bin/bash\n"
                for dir in dirs:
                    text += "cd %s\n"%pjoin(cwd,dir)
                    text += "./%s 1>./%s 2>&1 &\n"\
                        %(run_shower_dict[dir],"shower_screen_output.txt")
                text += "cd %s"%current_cwd
                offset = 0
                while offset<100000:
                    jobpath = pjoin(bin_dir,"showerjob%d"%(i+offset))
                    if not os.path.isfile(jobpath):
                        break
                    else:
                        offset += 1
                jobpath = pjoin(bin_dir,"showerjob%d"%(i+offset))
                jobfile = open(jobpath,"w")
                jobfile.write(text)
                jobfile.close()
                self.job_core["showerjob%d"%(i+offset)] = len(dirs)

            job_dict["bin"]= \
                [ file for file in os.listdir(bin_dir)\
                      if file.startswith('showerjob') ]

        ho_status = ['Showering events...' ]
        devnull = os.open(os.devnull, os.O_RDWR)
        if mode in ['HO_Shower']:
            for i,status in enumerate(ho_status):
                nevents_unweighted = []
                nevt_job = 0
                split = i == 1 and \
                    nevt_job >0 and \
                    any([int(l.split()[1])> int(nevt_job) \
                             for l in nevents_unweighted if l])
                self.update_status(status,level='shower')
                self.run_all(job_dict, [['2','B','%d' %i ]],status,split_jobs = split,cwd=cwd)

        if int(self.cluster_mode) == 1:
            # if cluster run, wait 10 sec so that event files are  transferred back
            self.update_status(\
                "Waiting while files are transferred back from the cluster nodes",\
                    level='shower')
            time.sleep(10)

        # combine generated topdrawer, gnuplot or root files
        self.collect_plots(cwd=cwd,dirs=run_shower_dirs,appstr="_"+shower_tag)

class ExtLauncher(object):
    """ Generic Class for executing external program """
    
    program_dir = ''
    executable = ''  # path from program_dir
    
    force = False

    def __init__(self, cmd, running_dir, card_dir='', **options):
        """ initialize an object """

        self.running_dir = running_dir
        self.card_dir = os.path.join(self.running_dir, card_dir)
        self.cmd_int = cmd
        #include/overwrite options 
        for key,value in options.items():
            setattr(self, key, value)
            
        self.cards = [] # files can be modified (path from self.card_dir) 

    def run(self):
        """ execute the main code """

        self.prepare_run()
        for card in self.cards:
            self.treat_input_file(card, default = 'n')

        self.launch_program()

    def prepare_run(self):
        """ aditional way to prepare the run"""
        pass

    def launch_program(self):
        """launch the main program"""
        subprocess.call([self.executable], cwd=self.running_dir)

    def edit_file(self, path):
        """edit a file"""

        path = os.path.realpath(path)
        misc.open_file(path)

    # Treat Nicely the timeout
    def timeout_fct(self,timeout):
        if timeout:
            # avoid to always wait a given time for the next answer
            self.force = True

    def ask(self, question, default, choices=[], path_msg=None):
        """nice handling of question"""

        if not self.force:
            return self.cmd_int.ask(question, default, choices=choices,
                                path_msg=path_msg, fct_timeout=self.timeout_fct)
        else:
            return str(default)

    def treat_input_file(self, filename, default=None, msg=''):
        """ask to edit a file"""

        if msg == '' and filename == 'param_card.dat':
            msg = \
            "WARNING: If you edit this file don\'t forget to consistently "+\
            "modify the different parameters,\n especially the width of all "+\
            "particles."

        if not self.force:
            if msg:  print msg
            question = 'Do you want to edit file: %(card)s?' % {'card':filename}
            choices = ['y', 'n']
            path_info = 'path of the new %(card)s' % {'card':os.path.basename(filename)}
            ans = self.ask(question, default, choices, path_info)
        else:
            ans = default

        if ans == 'y':
            path = os.path.join(self.card_dir, filename)
            self.edit_file(path)
        elif ans == 'n':
            return
        else:
            path = os.path.join(self.card_dir, filename)
            files.cp(ans, path)

class Pythia8Launcher(ExtLauncher):
    """A class to launch Pythia8 run"""

    def __init__(self, running_dir, cmd_int, **option):
        """ initialize launching Pythia 8"""
        
        running_dir = os.path.join(running_dir, 'examples')
        ExtLauncher.__init__(self, cmd_int, running_dir, '.', **option)
        self.cards = []

    def prepare_run(self):
        """ ask for pythia-pgs/delphes run """
        
        # Find all main_model_process.cc files
        date_file_list = []
        for file in glob.glob(os.path.join(self.running_dir,'main_*_*.cc')):
            # retrieves the stats for the current file as a tuple  
            # (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime)
            # the tuple element mtime at index 8 is the last-modified-date 
            stats = os.stat(file)
            # create tuple (year yyyy, month(1-12), day(1-31), hour(0-23), minute(0-59), second(0-59), 
            # weekday(0-6, 0 is monday), Julian day(1-366), daylight flag(-1,0 or 1)) from seconds since epoch
            # note:  this tuple can be sorted properly by date and time  
            lastmod_date = time.localtime(stats[8])
            date_file_list.append((lastmod_date, os.path.split(file)[-1]))

        if not date_file_list:
            raise HELACOniaError, 'No Pythia output found'

        # Sort files according to date with newest first
        date_file_list.sort()
        date_file_list.reverse()
        files = [d[1] for d in date_file_list]

        answer = ''
        answer = self.ask('Select a main file to run:', files[0], files)

        self.cards.append(answer)

        self.executable = self.cards[-1].replace(".cc","")

        # Assign a valid run name if not put in options
        if self.name == '':
            for i in range(1000):
                path = os.path.join(self.running_dir, '',
                                    '%s_%02i.log' % (self.executable, i))
                if not os.path.exists(path):
                    self.name = '%s_%02i.log' % (self.executable, i)
                    break

        if self.name == '':
            raise HELACOniaError, 'too many runs in this directory'

        # Find all exported models
        models = glob.glob(os.path.join(self.running_dir,os.path.pardir,
                                        "Processes_*"))
        models = [os.path.split(m)[-1].replace("Processes_","") for m in models]

        # Extract model name from executable 
        models.sort(key=len)
        models.reverse()
        model_dir = ""
        for model in models:
            if self.executable.replace("main_", "").startswith(model):
                model_dir = "Processes_%s" % model
                break
        if model_dir:
            self.model = model
            self.model_dir = os.path.realpath(os.path.join(self.running_dir,
                                                           os.path.pardir,
                                                           model_dir))
            self.cards.append(os.path.join(self.model_dir,
                                           "param_card_%s.dat" % model))

    def launch_program(self):
        """launch the main program"""

        # Make pythia8
        print "Running make for pythia8 directory"
        misc.compile(cwd=os.path.join(self.running_dir, os.path.pardir), mode='cpp')
        if self.model_dir:
            print "Running make in %s" % self.model_dir
            misc.compile(cwd=self.model_dir, mode='cpp')
        # Finally run make for executable
        makefile = self.executable.replace("main_","Makefile_")
        print "Running make with %s" % makefile
        misc.compile(arg=['-f', makefile], cwd=self.running_dir, mode='cpp')

        print "Running " + self.executable
        
        output = open(os.path.join(self.running_dir, self.name), 'w')
        if not self.executable.startswith('./'):
            self.executable = os.path.join(".", self.executable)
        subprocess.call([self.executable], stdout = output, stderr = output,
                        cwd=self.running_dir)

        # Display the cross-section to the screen 
        path = os.path.join(self.running_dir, self.name)
        pydoc.pager(open(path).read())
        
        print "Output of the run is found at " + \
            os.path.realpath(os.path.join(self.running_dir, self.name))
        
        

if __name__ == '__main__':
    #print os.environ
    #print os.path.realpath(__file__)
    #print root_path
    #ho_dir="/Users/erdissshaw/Works/HELAC-Onia/test-cluster/test2"
    job1=HELACOniaRunCmd(None,{'pt_list':["15d0","25d0"]})
    #job1.run_exe("Helac-Onia",['0', 'born', '0'],'monitor',ho_dir)
    #time.sleep(1)
    #job1.wait_for_complete('monitor')
