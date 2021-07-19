import os
import re
import math 
import numpy as np
from scipy.interpolate import interp1d
from lhereader import LHEReader

import matplotlib.pyplot as plt
from matplotlib import ticker, cm

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
def main():
  # get access to the data
  cases = [
    ["spinzero", "Drell-Yan, "+r'$\beta$-independent | spin: 0', "spinzero_eff.pdf"],
    ["spinhalf", "Drell-Yan, "+r'$\beta$-independent | spin: 1/2', "spinhalf_eff.pdf"],
    ["spinone", "Drell-Yan, "+r'$\beta$-independent | spin: 1', "spinone_eff.pdf"],
    ["spinzero_beta", "Drell-Yan, "+r'$\beta$-dependent | spin: 0', "spinzero_beta_eff.pdf"],
    ["spinhalf_beta", "Drell-Yan, "+r'$\beta$-dependent | spin: 1/2', "spinhalf_beta_eff.pdf"],
    ["spinone_beta", "Drell-Yan, "+r'$\beta$-dependent | spin: 1', "spinone_beta_eff.pdf"],
    ["PF_spinzero", "Photon Fusion, "+r'$\beta$-independent | spin: 0', "PF_spinzero_eff.pdf"],
    ["PF_spinhalf", "Photon Fusion, "+r'$\beta$-independent | spin: 1/2', "PF_spinhalf_eff.pdf"],
    ["PF_spinone", "Photon Fusion, "+r'$\beta$-independent | spin: 1', "PF_spinone_eff.pdf"],
    ["PF_spinzero_beta", "Photon Fusion, "+r'$\beta$-dependent | spin: 0', "PF_spinzero_beta_eff.pdf"],
#    ["PF_spinhalf_beta", "Photon Fusion, "+r'$\beta$-dependent | spin: 1/2', "PF_spinhalf_beta_eff.pdf"],
    ["PF_spinone_beta", "Photon Fusion, "+r'$\beta$-dependent | spin: 1', "PF_spinone_beta_eff.pdf"] ]

  nruns = 70
  runs = [None]*nruns
  for r in list(range(0, nruns)):
    runs[r] = f"{r+1:02}"

  # prepare ntuples [mmm, gch, acceptance]
  gchs = []
  mmms = []
  effs = []
  
  bekin = np.linspace(0., 4000., 100)
  vekin = []
  beta = np.linspace(-5., 5., 50)
  veta =[]
  bpt = np.linspace(0., 4000., 100)
  vpt =[]
  bbeta = np.linspace(0., 1., 50)
  vbeta =[]
  
  # loop over all cases
  for case in cases:
    outdir = case[0]
    pdfname = case[2]
    
    for run in runs:
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
      mmms.append(mmm)
      gchs.append(gch)
            
      # get access to event file
      dfile = evdir+'/unweighted_events.lhe'
      os.system("gunzip "+dfile+".gz")
      reader = LHEReader(dfile)

      # loop over events
      vekin.clear()
      veta.clear()
      vpt.clear()
      vbeta.clear()
      nevtot = 0
      nevALICE = 0
      
      for iev, event in enumerate(reader):
        nevtot += 1
        
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
          
        if inALICE == True:
          nevALICE += 1
          
      # update efficiencies
      eff = float(nevALICE)/float(nevtot)
      effs.append(eff)
      print("Run ",run, "Mass ", mmm, "Charge ", gch, "Efficiency ", eff)

      # plot histograms
      fig, axs = plt.subplots(2,2)
      ttxt = "| "+case[1]+" | m: "+str(mmm)+" | g: "+str(gch)+" |"
      
      nekin, bekin, pekin = axs[0][0].hist(vekin, bins=bekin, histtype='step', density=True)
      axs[0][0].set_title(ttxt, loc='left')
      axs[0][0].set_xlabel("kinetic energy")
      axs[0][0].set_yscale("log")
      neta, beta, peta = axs[0][1].hist(veta, bins=beta, histtype='step', density=True)
      axs[0][1].set_xlabel("eta")
      axs[0][1].set_yscale("log")
      npt, bpt, ppt = axs[1][0].hist(vpt, bins=bpt, histtype='step', density=True)
      axs[1][0].set_xlabel("pt")
      axs[1][0].set_yscale("log")
      nbeta, bbeta, pbeta = axs[1][1].hist(vbeta, bins=bbeta, histtype='step', density=True)
      axs[1][1].set_xlabel("beta")
      axs[1][1].set_yscale("log")

      plt.tight_layout()
      plt.savefig(outdir+"/mm_"+run+".pdf")
      plt.close()

      # clean up
      os.system("gzip "+dfile)
      
    # create efficiency map
    x = np.array(gchs)
    y = np.array(mmms)
    z = np.array(effs)
    
    # reshape arrays
    ny = len(np.where(y==y[0])[0])
    nx = int(len(y)/ny)
    x = x.reshape((nx, ny))
    y = y.reshape((nx, ny))
    z = z.reshape((nx, ny))
    
    # plot efficiency map
    fig, ax = plt.subplots()
    ttxt = "| "+case[1]+" |"
    ax.set_title(ttxt)
    ax.set_xlabel("magnetic charge")
    ax.set_ylabel("mass")
    cs = ax.contourf(x, y, z, levels=[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0], locator=ticker.LinearLocator(), cmap=cm.PuBu_r)
    cbar = fig.colorbar(cs, pad=0.02)
    cbar.set_label('ALICE TPC acceptance', rotation=270, labelpad=6.)
    
    plt.savefig(outdir+"/efficiency.pdf")
    plt.close()

# -----------------------------------------------------------------------------
if __name__ == "__main__":

  main()

# -----------------------------------------------------------------------------
