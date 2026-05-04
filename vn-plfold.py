#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# ==============================================================================
def print_history():
    print('2026-May-04: Original version.')

# ==============================================================================

class InteractiveFoldPlot:
    def __init__(self, w, trs, trsc, fac, istr, cmap='gist_stern', hic=99., loc=1.0):
        self.w = w
        self.trs = trs
        self.trsc = trsc
        self.fac = fac
        self.istr = istr
        self.hic = hic
        self.loc = loc

        # Setup Figure with shared X axis
        self.fig, self.axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
        plt.subplots_adjust(hspace=0.08)

        # Extent defines wavelength range and 2nd phase repeat
        self.extent = [w[0], w[-1], 0, 2]

        # Plot 1: Standard Folded Spectra using 'gist_stern' colormap
        self.im1 = self.axes[0].imshow(self.trs, aspect='auto', extent=self.extent,
                                       origin='lower', cmap=cmap, interpolation='nearest')
        # Plot 2: Inverse of Standard Folded
        self.im2 = self.axes[1].imshow(self.trs, aspect='auto', extent=self.extent,
                                       origin='lower', cmap=cmap+'_r', interpolation='nearest')
        # Plot 3: Atmosphere/Phase Subtracted
        self.im3 = self.axes[2].imshow(self.trsc, aspect='auto', extent=self.extent,
                                       origin='lower', cmap=cmap, interpolation='nearest')

#        self.axes[0].set_title(f"Zooming/panning in the UPPER panel will result in Interactive Robust Scaling (1st/99th): {istr}")
        self.axes[0].set_title(f"The color scale automatically recalculates based on the visible data whenever you zoom or pan the UPPER graph.\
        \nCheck also the parameters -fac -cmap -hic and -loc for ignoring outliers and finer color scaling.")
        self.axes[0].set_ylabel("Phase")
        self.axes[1].set_ylabel("Phase (Inverse)")
        self.axes[2].set_ylabel("Phase (Subtracted)")
        self.axes[2].set_xlabel("Wavelength")

        # Connect to zoom/pan events
        self.axes[0].callbacks.connect('xlim_changed', self.on_zoom)
        self.update_normalization()
        plt.show()

    def update_normalization(self):
        """Independently calculates robust scaling for each plot using finite visible pixels."""
        xlim = self.axes[0].get_xlim()
        idx = np.where((self.w >= xlim[0]) & (self.w <= xlim[1]))[0]

        if len(idx) > 5:
            # Data slices for the visible wavelength window
            vis_trs = self.trs[:, idx]
            vis_trsc = self.trsc[:, idx]

            # Filter for finite numbers to avoid NaN/Inf errors
            finite_trs = vis_trs[np.isfinite(vis_trs)]
            finite_trsc = vis_trsc[np.isfinite(vis_trsc)]

            if finite_trs.size > 0 and finite_trsc.size > 0:
                # Robust scaling using 1st and 99th percentiles
                v1_min, v1_max = np.percentile(finite_trs, [self.loc, self.hic])
                v3_min, v3_max = np.percentile(finite_trsc, [self.loc, self.hic])

                # Apply contrast factor 'fac' independently
                self.im1.set_clim(v1_min, v1_min + (v1_max - v1_min) / self.fac)
                self.im2.set_clim(v1_min, v1_min + (v1_max - v1_min) / self.fac)
                self.im3.set_clim(v3_min, v3_min + (v3_max - v3_min) / self.fac)

            self.fig.canvas.draw_idle()

    def on_zoom(self, event_ax):
        self.update_normalization()

def process_and_plot(fcol=1, tcol=2, wcol=1, scol=2, dcol=0, cfile='specnames',
                     nbin=None, data_dir='', rebin=0, fac=1.5, cmap='gist_stern',
                     hic=99., loc=1.):
    """Processes spectra in-memory and launches the interactive viewer."""

    if not os.path.exists(cfile):
        print(f"File {cfile} not found.")
        return

    # 1. Load list of spectra
    with open(cfile, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    fnam, pha = [], []
    for line in lines:
        parts = line.split()
        fnam.append(os.path.join(data_dir, parts[fcol-1]))
        pha.append(float(parts[tcol-1]))

    pha = np.array(pha)

    # 2. Load and interpolate spectra
    spin = []
    w1 = None
    for i, name in enumerate(fnam):
        data = np.loadtxt(name)
        w, sp = data[:, wcol-1], data[:, scol-1]

        if i == 0:
            w1 = w
            sp_ref_mean = np.mean(sp)
            spin.append(sp)
        else:
            sp_interp = np.interp(w1, w, sp)
            spin.append(sp_interp / np.mean(sp_interp) * sp_ref_mean)

    spin = np.array(spin)

    # 3. Handle Phase Folding or Sorting
    if rebin:
        idx = np.argsort(pha % 1.0)
        trsp = spin[idx]
    else:
        if nbin is None: nbin = len(pha)
        trsp = np.zeros((nbin, len(w1)))
        counts = np.zeros(nbin)
        for i in range(len(pha)):
            idx = int((pha[i] % 1.0) * nbin)
            trsp[idx] += spin[i]
            counts[idx] += 1
        mask = counts > 0
        trsp[mask] = (trsp[mask].T / counts[mask]).T

    # 4. Preparation for Visualization (2-phase repeat)
    trs = np.tile(trsp, (2, 1))

    # SAFE Atmospheric subtraction
    nwo = len(w1)
    mid = slice(int(0.4*nwo), int(0.6*nwo))
    clev = trs[:, mid].mean(axis=1)

    # Avoid division by zero: replace 0 with a very small number
    clev_safe = np.where(clev == 0, 1e-10, clev)
    mean_safe = clev_safe.mean()

    # Spectrum divided by phase-dependence
    trsc = trs / (clev_safe[:, np.newaxis] / mean_safe)

    # 5. Launch Viewer
    InteractiveFoldPlot(w1, trs, trsc, fac, cfile, cmap, hic, loc)

#########################################################################

def print_header():
    print ("")
    print ("***********************************************************************************")
    print ("**                                vn-plfold.py                                   **")
    print ("**              A little utility to plot spectra folded over phase               **")
    print ("**                                                                               **")
    print ("**                                2026-May-04                                    **")
    print ("**                             Vitaly  Neustroev                                 **")
    print ("**   (based on Henk Spruit's IDL code converted to Python with help of Gemini)   **")
    print ("***********************************************************************************")

def usage():
    #print_header()
    "Usage function"
    print ("Usage: %s [options] FileName" % sys.argv[0])
    print (" ")
    print ("FileName is a name of file containing list of spectra and having at least 2 columns:")
    print ("SpecNames and Phases. Each SpecName must have at least 2 columns: wavelengths and fluxes.")
    print ("The wavelength scale is taken from the 1st file, the next files are interpolated onto this scale.")
    print ("")
    print ("Options (the default values are given as examples):")
    print ("     -h    : Help")
    print ("     -H    : Summary of changes")
    print ("")
    print ("     -dir  : prefix to be added to file names from list    (example: -dir=./)")
    print ("     -fcol : column number of file name in list of spectra (example: -fcol=1)")
    print ("     -tcol : column number of phase value                  (example: -tcol=2)")
#    print ("     -dcol : column number of phase interval               (example: -dcol=0)")
    print ("     -wcol : column number of wavelength in spectrum file  (example: -wcol=1)")
    print ("     -scol : column number of spectrum                     (example: -scol=2)")
    print ("")
#    print ("     -dph  : a fixed phase interval can be specified here for the exposure times (example: -dph=0)")
    print ("     -rebin: if rebin=1, rebin data into equidistant phase bins (example: -rebin=0)")
    print ("     -nbin : number of phase bins to be used. Default: number of spectra")
    print ("                                                           (example: -nbin=None)")
    print ("     -cmap : color map for image                           (example: -cmap=gist_stern)")
    print ("     -fac  : contrast factor for image                     (example: -fac=1.0)")
    print ("     -hic and -loc: high and low cut of intensity (percentiles) (example: -hic=99 -loc=1)")

#    print ("")
#    sys.exit(-1)


##########################################################################


if __name__ == '__main__':

    if len(sys.argv) == 1:
        print_header()
        usage()

    data_dir=''
    fcol=1
    tcol=2
    dcol=0
    wcol=1
    scol=2
    cfile=''
    nbin=None
    rebin=0
    fac=1.0
    hic=99.0
    loc=1.0
    cmap='gist_stern'

    for i in range(len(sys.argv)-1):
        CmdLinePar = sys.argv[i+1]
        if CmdLinePar[0] == '-':
            if CmdLinePar == '-h':
            #if ('h') in CmdLinePar[1:2]:
                usage()
                exit()
            if ('H') in CmdLinePar[1:2]:
                print_header()
                print ("\nUsage: %s [options] FileName" % sys.argv[0])
                print_history()
                exit()
            if ('dir=') in CmdLinePar[1:5]:
                data_dir = CmdLinePar[5:]
            if ('fcol=') in CmdLinePar[1:6]:
                fcol=int(CmdLinePar[6:])
            if ('tcol=') in CmdLinePar[1:6]:
                tcol=int(CmdLinePar[6:])
            if ('wcol=') in CmdLinePar[1:6]:
                wcol=int(CmdLinePar[6:])
            if ('scol=') in CmdLinePar[1:6]:
                scol=int(CmdLinePar[6:])
            if ('dcol=') in CmdLinePar[1:6]:
                dcol=int(CmdLinePar[6:])
            if ('dph=') in CmdLinePar[1:5]:
                dph=float(CmdLinePar[5:])
            if ('nbin=') in CmdLinePar[1:6]:
                nbin=int(CmdLinePar[6:])
            if ('rebin=') in CmdLinePar[1:7]:
                rbin=int(CmdLinePar[7:])
                if rbin != 0: rebin=1
                else:rebin=0
            if ('cmap=') in CmdLinePar[1:6]:
                cmap=CmdLinePar[6:]
            if ('fac=') in CmdLinePar[1:5]:
                fac=float(CmdLinePar[5:])
            if ('hic=') in CmdLinePar[1:5]:
                hic=float(CmdLinePar[5:])
            if ('loc=') in CmdLinePar[1:5]:
                loc=float(CmdLinePar[5:])
        else:
            cfile = CmdLinePar
            #print("Name: ",cfile)
            if not os.path.isfile(cfile):
                print("The File ",cfile," doesn't exist.")
                #exit(-1)
    if cfile=='':
        cfile = input("\nEnter the file name of a list of spectra: ")
        if cfile=='': exit(0)
        elif not os.path.isfile(cfile):
            print("The File ",cfile," doesn't exist.")
            exit(-1)
##            FileName = input("Enter the file name of a list of spectra: ")
##        try:
##            infile = open(FileName[s1:], "r")
##            lines = infile.readlines()
##            infile.close()
##            print("The FileList ",FileName[s1:]," includes ",len(lines)," FileNames")
##            FileList=True
##        except:
##            print("Something wrong with the FileList ",FileName[s1:])

    process_and_plot(fcol=fcol, tcol=tcol, wcol=wcol, scol=scol, dcol=dcol, cfile=cfile,
        data_dir=data_dir, rebin=rebin, nbin=nbin, fac=fac*0.7, cmap=cmap, hic=hic, loc=loc)
