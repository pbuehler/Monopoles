import os
import re
import math 
import numpy as np
from scipy.interpolate import interp1d
from lhereader import LHEReader

import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import matplotlib.backends.backend_pdf

# global definition of minimum kinetic energy necessary for Monopole to reach TPC
gDs = [0, 1, 2, 3, 4, 5, 6, 7]
Emins = [0., 20., 50., 90., 150., 220., 320., 450. ]
iEkinmin = interp1d(gDs, Emins)

# -----------------------------------------------------------------------------
def isInALICE(mm, gch):
  if mm.p4().eta > -0.9 and  mm.p4().eta < 0.9:
    ekin = math.sqrt(mm.p4().m*mm.p4().m+mm.p4().p*mm.p4().p)-mm.p4().m
    if ekin > ekinminALICE(gch):
      return True
    else:
      return False
  else:
    return False

# -----------------------------------------------------------------------------
def ekinminALICE(gch):
  return iEkinmin(gch)

# -----------------------------------------------------------------------------
def main(title, cases, pdfname):
  # prepare data buffers
  bekin = np.linspace(0., 4000., 100)
  vekin = []
  beta = np.linspace(-5., 5., 50)
  veta =[]
  bpt = np.linspace(0., 4000., 100)
  vpt =[]
  bbeta = np.linspace(0., 1., 50)
  vbeta =[]
  
  # initialize plots
  fig1, ax1 = plt.subplots()
  ax1.set_title(title)
  ax1.set_xlabel("kinetic energy")
  ax1.set_yscale("log")
  ax1.set_xlim([0., 4000.])
  
  fig2, ax2 = plt.subplots()
  ax2.set_title(title)
  ax2.set_xlabel("eta")
  ax2.set_yscale("log")
  ax2.set_xlim([-5., 5.])

  fig3, ax3 = plt.subplots()
  ax3.set_title(title)
  ax3.set_xlabel("pt")
  ax3.set_yscale("log")
  ax3.set_xlim([0., 4000.])

  fig4, ax4 = plt.subplots()
  ax4.set_title(title)
  ax4.set_xlabel("beta")
  ax4.set_yscale("log")
  ax4.set_xlim([0., 1.])
  
  # loop over all cases
  for case in cases:
    outdir = case[0]
    run = case[1]
    evdir = outdir+'/Events/run_'+run
    
    # get run parameters from banner file
    #   mass
    #   magnetic charge
    bfname = evdir+"/run_"+run+"_tag_1_banner.txt"
    mmm = -1.
    gch = -1.
    f = open(bfname, "r")
    mmmfound = False
    gchfound = False
    for line in f:
      if re.search("# mmm",line):
        mmm = float(line.split()[1])
      if re.search("# gch",line):
        gch = float(line.split()[1])
    f.close()
    
    if (mmm < 0.) | (gch < 0.):
      raise Exception("Couldn't properly parse %s", bfname)
            
    # get access to event file
    dfile = evdir+'/unweighted_events.lhe'
    os.system("gunzip "+dfile+".gz")
    reader = LHEReader(dfile)

    # loop over events
    vekin.clear()
    veta.clear()
    vpt.clear()
    vbeta.clear()
      
    for iev, event in enumerate(reader):      
      # Find MM particles
      mms = filter(lambda x: abs(x.pdgid)==4110000, event.particles)
      inALICE = False
      for mm in mms:    
        inALICE |= isInALICE(mm, gch)
        
        ekin = math.sqrt(mm.p4().m*mm.p4().m+mm.p4().p*mm.p4().p)-mm.p4().m
        vekin.append(ekin)
        veta.append(mm.p4().eta)
        vpt.append(mm.p4().pt)
        vbeta.append(mm.p4().beta)
        
    print("Outdir ", outdir, "Run ",run, "Mass ", mmm, "Charge ", gch)
    os.system("gzip "+dfile)

    # plot histograms
    n, b, p = ax1.hist(vekin, bins=bekin, histtype='step', density=True, label=case[2])
    n, b, p = ax2.hist(veta, bins=beta, histtype='step', density=True, label=case[2])
    n, b, p = ax3.hist(vpt, bins=bpt, histtype='step', density=True, label=case[2])
    n, b, p = ax4.hist(vbeta, bins=bbeta, histtype='step', density=True, label=case[2])

  legend = ax1.legend(loc='upper right', shadow=True, fontsize='small')
  legend = ax2.legend(loc='lower center', shadow=True, fontsize='small')
  legend = ax3.legend(loc='upper right', shadow=True, fontsize='small')
  legend = ax4.legend(loc='upper left', shadow=True, fontsize='small')

  pdf = matplotlib.backends.backend_pdf.PdfPages(pdfname)
  pdf.savefig( fig1 )
  pdf.savefig( fig2 )
  pdf.savefig( fig3 )
  pdf.savefig( fig4 )
  pdf.close()

  plt.close()

# -----------------------------------------------------------------------------
if __name__ == "__main__":

  title = "Photon Fusion, "+r'$\beta$-independent'
  cases = [
    ["PF_spinzero", "15", "spin 0, m=1000, 1gD"],
    ["PF_spinhalf", "15", "spin 1/2, m=1000, 1gD"],
    ["PF_spinone",  "15", "spin 1, m=1000, 1gD"] ]
  pdfname = "PF_1_1000_beta.pdf"
    
  main(title, cases, pdfname)

# -----------------------------------------------------------------------------
