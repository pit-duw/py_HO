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
import glob
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

combine_factors = []
gener = 1 # 0 phegas 1 vegas
pt_list = []

def collect_results(cwd=None):
    """collect the results for parallel computations."""
    p_dirs = [file for file in os.listdir(cwd)
              if file.startswith('P') and \
                  os.path.isdir(pjoin(cwd,file))]
    result_dir = pjoin(cwd,"results")
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    ptQ=False
    if len(pt_list) == 0:
        # total cross section
        text = "# Total Cross Section"
        text += "\n# SubDir            NMAX    ACC                       SD                      "
    else:
        # pt distribution
        text = "# Pt Differential Distribution"
        text += "\n# SubDir            NMAX    Pt          ACC                       SD                      "
        ptQ = True

    for i,p_dir in enumerate(p_dirs):
        tmp_dir = p_dir.strip()
        tmp_len = len(tmp_dir)
        np_dir = " "*max(17,tmp_len)
        np_dir = tmp_dir+np_dir[tmp_len:]
        screenoutput=pjoin(cwd,p_dir,"output","screen_output.txt")
        tmpnmax, tmpacc, tmpsd =read_result(screenoutput,gener=gener)
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
            tmppt=pt_list[i].lower().replace("d","e").strip()
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
    return result_path

def collect_plots(cwd=None):
    """collect the plots for parallel computations."""
    p_dirs = [file for file in os.listdir(cwd)
              if file.startswith('P') and \
                  os.path.isdir(pjoin(cwd,file))]
    result_dir = pjoin(cwd,"results")
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    if not p_dirs:return
    td_plots = []
    gnu_plots = []
    root_plots = []
    for i,p_dir in enumerate(p_dirs):
        td_plots.extend([pjoin(cwd,p_dir,"output",file) \
                             for file in os.listdir(pjoin(cwd,p_dir,"output"))\
                             if file.endswith('.top') and \
                             os.path.isfile(pjoin(cwd,p_dir,"output",file))])
        gnu_plots.extend([pjoin(cwd,p_dir,"output",file) \
                              for file in os.listdir(pjoin(cwd,p_dir,"output"))\
                              if file.endswith('.gnu') and \
                              os.path.isfile(pjoin(cwd,p_dir,"output",file))])
        root_plots.extend([pjoin(cwd,p_dir,"output",file) \
                               for file in os.listdir(pjoin(cwd,p_dir,"output"))\
                               if file.endswith('.C') and \
                               os.path.isfile(pjoin(cwd,p_dir,"output",file))])

    # combine topdrawer file
    if td_plots:
        combine_td_gnu(td_plots,result_dir,type="td")
    if gnu_plots:
        combine_td_gnu(gnu_plots,result_dir,type="gnu")
    if root_plots:
        combine_root(root_plots,result_dir)

def combine_root(root_plots,result_dir):
    """combine root files"""
    rt_text = ""
    rt_combine_factors = combine_factors
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
                if plot_start(line,type="root"):
                    store_res = 1
                    i = i + 1
                elif plot_end(line,type="root"):
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
                elif plot_start(line,type="root"):
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
                        if dim_type == 2:
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
                if plot_start(line,type="root"):
                    store_res = 1
                    i = i + 1
                elif plot_end(line,type="root"):
                    store_res = 0

                if store_res == 1 and not plot_start(line,type="root"):
                    if 'xybin' in line:
                        dim_type=3
                    else:
                        dim_type=2
                    rep_str = 'plot_%i'%i
                    data = line.split()
                    data.pop()
                    if 'FindBin' in line:
                        if dim_type==2:
                            j=find_pos(float(data[-1]),store_plots[rep_str])
                        else:
                            j=find_pos_3d([float(data[-3]),float(data[-1])],store_plots[rep_str])
                    if 'SetBinContent' in line:
                        store_plots[rep_str][j][1]+=float(data[-1])*yfactor
                    if 'SetBinError' in line:
                        store_plots[rep_str][j][2]+=float(data[-1])*yfactor

        rt_file.close()

    for key,value in store_plots.items():
        if dim_type_list[store_index[key]-1]==2:
            store_plots_str[key]=root_print(value,i=store_index[key])
        else:
            store_plots_str[key]=root_print_3d(value,i=store_index[key])

    rt_text = rt_text%store_plots_str
    res_filename = "results.C"
    result_path = pjoin(result_dir,res_filename)
    result_file = open(result_path,"w")
    result_file.write(rt_text)
    result_file.close()                

def find_pos(x0,value_list):
    j=0
    for x,y,dy in value_list:
        if abs(x0-x)/max(abs(x0),abs(x),1e-99) < 1e-5:
            return j
        if x0 < x: break
        j=j+1
    value_list.insert(j,[x0,0.,0.])

    return j

def find_pos_3d(xy0,value_list):
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
    

def combine_td_gnu(td_plots,result_dir,type="td"):
    """combine topdrawer/gnuplot files"""
    td_text = ""
    td_combine_factors = combine_factors
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
                if plot_start(line,type=type):
                    store_res = 1
                    i = i + 1
                    dim_type_list.append(plot_dimtype(line,type=type))
                elif plot_end(line,type=type):
                    store_res = 0

                        
                if store_res == 0:
                    if '\n' in line:
                        td_text +='%s'%line
                    else:
                        td_text +='%s\n'%line
                elif plot_start(line,type=type):
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
                        elif type=='gnu':
                            data[2] = data[2]*yfactor
                            store_plots[rep_str].append(data)
                    elif dim_type_list[i-1] == 3 and type == "gnu":
                        store_plots[rep_str].append(None)

        else:
            j = 0
            for line in td_file:
                if plot_start(line,type=type):
                    store_res = 1
                    i = i + 1
                elif plot_end(line,type=type):
                    store_res = 0
                    j=0

                if store_res == 1 and not plot_start(line,type=type):
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
        store_plots_str[key]=plot_print(value,type=type)

    if type=="gnu":
        td_text=td_text.replace('set format y "10^{%T}"','set format y "10^{%%T}"')
    td_text = td_text%store_plots_str
    res_filename = "results"
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

def plot_dimtype(line,type="td"):
    if type == "td":
        # only 2-dim for topdrawer
        return 2
    elif type == "gnu":
        if 'splot "-" with pm3d' in line:
            return 3
        else:
            return 2

def plot_start(line,type="td"):
    if type == "td":
        return td_start(line)
    elif type == "gnu":
        return gnu_start(line)
    elif type == "root":
        return root_start(line)
    else:
        return False

def plot_end(line,type="td"):
    if type == "td":
        return td_end(line)
    elif type == "gnu":
        return gnu_end(line)
    elif type == "root":
        return root_end(line)
    else:
        return False

def root_start(line):
    if 'SetStats(false)' in line:
        return True
    else:
        return False

def root_end(line):
    if 'histos -> Add(hist);' in line or 'hist_3d -> SetOption("colz");' in line:
        return True
    else:
        return False

def gnu_start(line):
    if 'plot "-" with histeps' in line or 'splot "-" with pm3d' in line:
        return True
    else:
        return False

def gnu_end(line):
    newline=line.split()
    if len(newline) == 1 and newline[0] == 'e':
        return True
    else:
        return False

def td_start(line):
    if 'SET ORDER X Y DY' in line:
        return True
    else:
        return False

def td_end(line):
    if 'HISTO' in line:
        return True
    else:
        return False

def plot_print(value_list,type="td"):
    """output topdrawer/gnuplot file"""
    if type == "td":
        return td_print(value_list)
    elif type == "gnu":
        return gnu_print(value_list)
    
    return ""

def td_print(value_list):
    """output topdrawer file"""
    text = ""
    for data in value_list:
        text += "       %5.4e       %5.4e       %5.4e\n"%tuple(data)
        
    return text

def gnu_print(value_list):
    """output gnuplot file"""
    text=""
    for data in value_list:
        if data==None:
            text +="     \n"
        else:
            text +="     %7.6e       %7.6e    %7.6e\n"%tuple(data)

    return text

def root_print(value_list,i=1):
    """output root file"""
    text=""
    for data in value_list:
        text +='   int xbin = id%i->FindBin(   %17.16e      );\n'%(i,data[0])
        text +='  id%i -> SetBinContent( xbin,   %17.16e );\n'%(i,data[1])
        text +='  id%i -> SetBinError( xbin,    %17.16e      );\n'%(i,data[2])

    return text

def root_print_3d(value_list,i=1):
    """output root file"""
    text=""
    for data in value_list:
        text +='   int xybin = id3d%i->FindBin(   %17.16e      ,   %17.16e      );\n'%(i,data[0][0],data[0][1])
        text +='  id3d%i -> SetBinContent( xybin,   %17.16e );\n'%(i,data[1])
        text +='  id3d%i -> SetBinError( xybin,    %17.16e      );\n'%(i,data[2])

    return text

def read_result(file,gener=1):
    """read the result from file.
    gener:0 phegas; 1, vegas."""
    if not isinstance(file,str) or not os.path.isfile(file):
        raise TypeError("%s is not a valid file"%str(file))
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
            # PHEGAS
            pass
        else:
            raise TypeError("unknown type of gener = %d"%gener)
        preline = line
    result_file.close()
    res_list.reverse()
    res_nmax = "nan"
    res_acc = "nan"
    res_sd = "nan"
    if len(res_list) == 0:
        print "no result in %s"%file
    for ncal,acc,sd in res_list:
        res_nmax = ncal
        res_acc = acc
        res_sd = sd
        if not "nan" in acc.lower() and not "nan" in sd.lower():
            break
    else:
        print "no valid result in %s"%file
    return res_nmax,res_acc,res_sd

def convertgnu2root(gnuplot,rootfile):
    """convert gnuplot file into root file"""
    rt_text = \
"""
 {
 hohisto = new TFile("results.root","recreate");
 hohisto -> cd();
 histos = new TObjArray(0);
"""
    store_plots_str = {}
    store_plots = {}
    store_index = {}

    rt_file = open(gnuplot,'r')
    store_res = 0
    i = 0
    title_pat = re.compile(r'\s*set title\s*\"(?P<title>[\s\S]*)\"\s*font\s*[\s\S]*')
    xlabel_pat = re.compile(r'\s*set xlabel\s*\"(?P<xlabel>[\s\S]*)\"\s*font\s*[\s\S]*')
    ylabel_pat = re.compile(r'\s*set ylabel\s*\"(?P<ylabel>[\s\S]*)\"\s*font\s*[\s\S]*')
    zlabel_pat = re.compile(r'\s*set zlabel\s*\"(?P<zlabel>[\s\S]*)\"\s*font\s*[\s\S]*')
    xrange_pat = re.compile(r'\s*set xrange\s*\[\s*(?P<xlow>[-]*\d*.\d*)\s*:\s*(?P<xup>[-]*\d*.\d*)\s*\][\s\S]*')
    yrange_pat = re.compile(r'\s*set yrange\s*\[\s*(?P<ylow>[-]*\d*.\d*)\s*:\s*(?P<yup>[-]*\d*.\d*)\s*\][\s\S]*')
    replace = {}
    threedplot=False
    nthreedplot=0
    threedplotdata = []
    threedplotept = False
    threedycount = 0
    threeddy = 0.
    threeddx = 0.
    threedxcount = 0
    for line in rt_file:
        if "set title" in line:
            title_match = title_pat.match(line)
            title = title_match.group('title')
            replace['title']=title
        elif "set xlabel" in line:
            xlabel_match = xlabel_pat.match(line)
            xlabel = xlabel_match.group('xlabel')
            replace['xlabel']=xlabel
        elif "set ylabel" in line:
            ylabel_match = ylabel_pat.match(line)
            ylabel = ylabel_match.group('ylabel')
            replace['ylabel']=ylabel
            nthreedplot=nthreedplot+1
        elif "set zlabel" in line:
            zlabel_match = zlabel_pat.match(line)
            zlabel = zlabel_match.group('zlabel')
            replace['zlabel']=zlabel
        elif "set xrange" in line:
            xrange_match = xrange_pat.match(line)
            replace['xlow']=xrange_match.group('xlow')
            replace['xup']=xrange_match.group('xup')
        elif "set yrange" in line:
            yrange_match = yrange_pat.match(line)
            replace['ylow']=yrange_match.group('ylow')
            replace['yup']=yrange_match.group('yup')
            threedplot=True
        elif plot_start(line,type='gnu'):
            store_res = 1
            i = i + 1
            if not threedplot:
                nxbin = 0
                replace['id']='id%i'%i
                rt_text_temp =\
"""
 hohisto -> cd();
 TH1F *hist = new TH1F( "%(id)s","%(title)s",%(nxbin)i,%(xlow)s,%(xup)s);

 %(id)s -> GetXaxis() -> SetTitle("%(xlabel)s");
 %(id)s -> GetYaxis() -> SetTitle("%(ylabel)s");

 %(id)s -> GetYaxis() -> SetTitleOffset(1.2);
 %(id)s -> SetStats(false);
"""
            else:
                nxbin = 0
                nybin = 0
                threedplotdata = []
                threedycount = 0
                threeddy = 0.
                threedxcount = 0
                threeddx = 0.
                replace['id']='id3d%i'%i
                if nthreedplot == 1:
                    rt_text_temp =\
"""

 Int_t colors[50];
 Int_t number = 3;
 Double_t red[number] = { 1.00, 0.00, 0.00 };
 Double_t green[number] = { 0.00, 1.00, 0.00 };
 Double_t blue[number] = { 1.00, 0.00, 1.00 };
 Double_t length[number] = { 0.00, 0.50, 1.00 };
 Int_t nb=50;
 Int_t fi = TColor::CreateGradientColorTable(number,

length,red,green,blue,nb);
 for (int i=0;i<50;i++) colors[i]=fi+i;
 gStyle -> SetPalette(50,colors);

"""
                else:
                    rt_text_temp = ""
                rt_text_temp += \
"""
 hohisto -> cd();
 TH2F *hist_3d = new TH2F( "%(id)s", "%(title)s", %(nxbin)i, %(xlow)s, %(xup)s, %(nybin)i, %(ylow)s, %(yup)s );

 %(id)s -> GetXaxis() -> SetTitle("%(xlabel)s");
 %(id)s -> GetYaxis() -> SetTitle("%(ylabel)s");
 %(id)s -> GetZaxis() -> SetTitle("%(zlabel)s");

 %(id)s -> GetZaxis() -> SetTitleOffset(1.2);
 %(id)s -> SetStats(false);
 
"""
        elif plot_end(line,type='gnu'):
            store_res = 0
            if not threedplot:
                rt_text_temp +="\n histos -> Add(hist);\n"
            else:
                rt_text_temp +="\n hist_3d -> SetOption(\"colz\");"
                rt_text_temp +="\n histos -> Add(hist_3d);\n"
                replace['nybin']=nybin
                nxbin = threedxcount/2
            replace['nxbin']=nxbin
            rt_text_temp = rt_text_temp%replace
            threedplot=False
            threedplotept = False
            rt_text += "\n\n"+rt_text_temp

        if store_res == 1 and not plot_start(line,type='gnu'):
            if not threedplot:
                data = line.split()
                nxbin = nxbin + 1
                if float(data[1]) > 0:
                    rt_text_temp +="\n int xbin = %(id)s->FindBin( "+data[0]+" );"
                    rt_text_temp +="\n %(id)s -> SetBinContent( xbin, "+data[1]+" );"
                    rt_text_temp +="\n %(id)s -> SetBinError( xbin, "+data[2]+" );"
            else:
                data = line.split()
                if not threedplotept and data:
                    threedycount = threedycount+1
                    if threedycount == 2:
                        threeddy = float(data[1])-float(threedplotdata[0][1])
                        threeddy = threeddy/2.0
                    threedplotdata.append(data)
                elif not data and not threedplotept:
                    nybin=threedycount/2
                    threedplotept=True
                    threedycount=0
                elif threedplotept and threedxcount == 0 and data:
                    threedxcount=2
                    threedycount=threedycount+1
                    threeddx= float(data[0])-float(threedplotdata[0][0])
                    threeddx=threeddx/2.0
                    rt_text_temp +="\n int xybin = %(id)s->FindBin( "+\
                        ("%17.16e"%(float(data[0])-threeddx))+" , "+\
                        ("%17.16e"%(float(data[1])+threeddy))+" );"
                    rt_text_temp +="\n %(id)s -> SetBinContent( xybin, "+data[2]+" );"
                    rt_text_temp +="\n %(id)s -> SetBinError(xybin, 0.0000000000000000e+00      );"
                elif threedplotept and data:
                    threedycount=threedycount+1
                    if threedycount%2 == 1 and threedxcount%2 == 0 and float(data[2]) > 0:
                        rt_text_temp +="\n int xybin = %(id)s->FindBin( "+\
                            ("%17.16e"%(float(data[0])-threeddx))+" , "+\
                            ("%17.16e"%(float(data[1])+threeddy))+" );"
                        rt_text_temp +="\n %(id)s -> SetBinContent( xybin, "+data[2]+" );"
                        rt_text_temp +="\n %(id)s -> SetBinError(xybin, 0.0000000000000000e+00      );"
                elif threedplotept and not data:
                    threedycount = 0
                    threedxcount = threedxcount + 1
    rt_text +=\
"""

 hohisto -> cd();
 if (histos -> GetEntries() > 0) then {
   histos->Write();
   hohisto -> Close();
 }
}

"""
    rt_file.close()
    result_file = open(rootfile,'w')
    result_file.write(rt_text)
    result_file.close()
    return
                    
if __name__ == '__main__':
    combine_factors = [1.,1.,2.]
    collect_plots(cwd="/Users/erdissshaw/Works/Plots/JpsiJpsi/ATLAS_20150227/PROC_HO_5")
    #convertgnu2root("/Users/erdissshaw/Works/Plots/JpsiJpsi/ATLAS_20150227/PROC_HO_1/results/results.gnu","results.C")
    #print os.environ
    #print os.path.realpath(__file__)
    #print root_path
    #ho_dir="/Users/erdissshaw/Works/HELAC-Onia/test-cluster/test2"
    #job1=HELACOniaRunCmd(None,{'pt_list':["15d0","25d0"]})
    #job1.run_exe("Helac-Onia",['0', 'born', '0'],'monitor',ho_dir)
    #time.sleep(1)
    #job1.wait_for_complete('monitor')
