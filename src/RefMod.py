# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 16:55:52 2022

@author: alaguillog
"""

import warnings
warnings.filterwarnings(
    "ignore",
    message="Warning: OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning
)
from ast import literal_eval
import argparse
from bisect import bisect_left
import concurrent.futures
import configparser
from datetime import datetime
import glob
import itertools
import logging
import math
import numpy as np
import os
import pandas as pd
from pathlib import Path
import pyopenms
import re
import sys
import shutup
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'
shutup.please()

def checkParams(mass, infiles):
    min_frag_mz = int(mass._sections['Spectrum Processing']['min_fragment_mz'])
    max_frag_mz = int(mass._sections['Spectrum Processing']['max_fragment_mz'])
    if (max_frag_mz > 0) & (max_frag_mz <= min_frag_mz):
        logging.error('max_frag_mz must be either 0 or a value greater than min_frag_mz')
        return(1)
    prot_column = str(mass._sections['FDR']['prot_column'])
    for f in infiles:
        cols = pd.read_csv(f, index_col=0, nrows=0, sep="\t").columns.tolist()
        if prot_column not in cols:
            logging.error('The file ' + str(f) + ' does not contain a ' + str(prot_column) + ' column. Please check the name of the protein column.')
            return(1)
    return(0)

def preProcess(args):
    if os.path.isdir(args.infile):
        infiles = []
        for f in os.listdir(args.infile):
            if f.lower().endswith(".tsv"):
                if len(args.dia) > 0:
                    infiles += [os.path.join(args.infile, os.path.basename(f).split(sep=".")[0] + "_ch" + str(c) + ".tsv") for c in args.dia]
                else:
                    infiles += [os.path.join(args.infile, f)]
        if len(infiles) == 0:
            sys.exit("ERROR: No TSV files found in directory " + str(args.infile))
    elif '*' in args.infile:
        infiles_t = glob.glob(args.infile)
        infiles = []
        if len(args.dia) > 0:
            for f in infiles_t:
                infiles += [os.path.join(args.infile, os.path.basename(f).split(sep=".")[0] + "_ch" + str(c) + ".tsv") for c in args.dia]
        else:
            infiles = infiles_t
    else:
        if len(args.dia) > 0:
            infiles = [os.path.join(os.path.dirname(args.infile), os.path.basename(args.infile).split(sep=".")[0] + "_ch" + str(c) + ".tsv") for c in args.dia]
        else:
            infiles = [args.infile]

    if os.path.isdir(args.rawfile):
        rawfiles = []
        rawbase = []
        for f in os.listdir(args.rawfile):
            if (f.lower().endswith(".mzml") or f.lower().endswith(".mgf")):
                rawfiles += [os.path.join(args.rawfile, f)]
                rawbase += [os.path.basename(f).split(sep=".")[0]]
        if len(rawfiles) == 0:
            sys.exit("ERROR: No mzML or MGF files found matching pattern " + str(args.rawfile))
    else:
        rawfiles = [args.rawfile]
        rawbase = [os.path.basename(args.rawfile).split(sep=".")[0]]

    return(infiles, rawfiles, rawbase)

def readRaw(msdata):
    if os.path.splitext(msdata)[1].lower() == ".mzml":
        mode = "mzml"
        logging.info("Reading mzML file (" + str(os.path.basename(Path(msdata))) + ")...")
        fr_ns = pyopenms.MSExperiment() # TODO: use OnDiscMSExperiment() to manage memory usage
        pyopenms.MzMLFile().load(str(msdata), fr_ns)
        index2 = 0
        logging.info("\t" + str(fr_ns.getNrSpectra()) + " spectra read.")
    elif os.path.splitext(msdata)[1].lower() == ".mgf":
        mode = "mgf"
        logging.info("Reading MGF file (" + str(os.path.basename(Path(msdata))) + ")...")
        fr_ns = pd.read_csv(msdata, header=None)
        index2 = fr_ns.to_numpy() == 'END IONS'
        logging.info("\t" + str(sum(fr_ns[0].str[:4]=="SCAN")) + " spectra read.")
        # logging.info("\t" + str(fr_ns[0].str.count('SCANS').sum()) + " spectra read.") # VERY SLOW
    else:
        logging.info("MS Data file type '" + str(os.path.splitext(msdata)[1]) + "' not recognized!")
        sys.exit()
    return(fr_ns, mode, index2)

def deisotope(ions, m_proton, max_charge):
    # Sort by descending intensity
    ions = np.array([ions[0][ions[1].argsort()][::-1], ions[1][ions[1].argsort()][::-1]])
    if len(ions[0]) > 1000: # Remove low intensity peaks
        ions = np.array([ions[0][:999], ions[1][:999]])
    # Sort by ascending m/z
    ions = np.array([ions[0][ions[0].argsort()], ions[1][ions[0].argsort()]])
    isoids = []
    for i in range(1, max_charge): # Consider charges up to precursor charge state - 1
        isoids1 = [np.isclose(ions[0], j-(m_proton/i), atol=0.005, rtol=0).any() for j in ions[0]]
        isoids2 = [np.isclose(ions[0], j-((m_proton*2)/i), atol=0.005, rtol=0).any() for j in ions[0]]
        isoids3 = [np.isclose(ions[0], j-((m_proton*2)/i), atol=0.005, rtol=0).any() for j in ions[0]]
        isoids += [np.logical_or.reduce((isoids1, isoids2, isoids3))]
    isoids = np.logical_or.reduce(isoids)
    ions = np.array([np.delete(ions[0], isoids), np.delete(ions[1], isoids)])
    return(ions)

def locateScan(scan, mode, fr_ns, spectra, spectra_n, index2, top_n, bin_top_n, min_ratio,
               min_frag_mz, max_frag_mz, m_proton, deiso, max_charge):
    if mode == "mgf":
        # index1 = fr_ns.to_numpy() == 'SCANS='+str(int(scan))
        try:
            index1 = fr_ns.loc[fr_ns[0]=='SCANS='+str(scan)].index[0] + 1
            # index1 = np.where(index1)[0][0]
        except IndexError:
            logging.info("\tERROR: Scan number " + str(scan) + " not found in MGF file.")
            sys.exit()
        index3 = np.where(index2)[0]
        index3 = index3[np.searchsorted(index3,[index1,],side='right')[0]]
        try:
            ions = fr_ns.iloc[index1:index3,:]
            ions[0] = ions[0].str.strip()
            ions[['MZ','INT']] = ions[0].str.split(" ",expand=True,)
            ions = ions.drop(ions.columns[0], axis=1)
            ions = ions.apply(pd.to_numeric)
        except ValueError:
            ions = fr_ns.iloc[index1+4:index3,:]
            ions[0] = ions[0].str.strip()
            ions[['MZ','INT']] = ions[0].str.split(" ",expand=True,)
            ions = ions.drop(ions.columns[0], axis=1)
            ions = ions.apply(pd.to_numeric)
        ions = np.array(ions.T)
        ions0 = ions[0]
        ions1 = ions[1]
    elif mode == "mzml":
        try:
            s = spectra[spectra_n.index(scan)]
        except AssertionError or OverflowError:
            logging.info("\tERROR: Scan number " + str(scan) + " not found in mzML file.")
            sys.exit()
        peaks = s.get_peaks()
        ions0 = peaks[0]
        ions1 = peaks[1]
    # Normalize intensity
    ions1 = (ions1/max(ions1))*100
    # Remove peaks below min_ratio
    if min_ratio > 0:
        cutoff1 = ions1/max(ions1) >= min_ratio
        ions0 = ions0[cutoff1]
        ions1 = ions1[cutoff1]
    # Return only top N peaks
    if bin_top_n:
        bins = np.digitize(ions[0], np.arange(55,max(ions[0]),110))
        ions_f = []
        for i in np.unique(bins):
            ions_t = np.array([ions[0][np.where(bins==i)], ions[1][np.where(bins==i)]])
            cutoff = len(ions_t[0])-top_n
            if (cutoff < 0) or (cutoff >= len(ions_t[0])): cutoff = 0
            ions_f += [np.array([ions_t[0][ions_t[1].argsort()][cutoff:], ions_t[1][ions_t[1].argsort()][cutoff:]])]
        ions = np.concatenate(ions_f, axis=1)
    elif top_n > 0:
        cutoff1 = ions1 >= ions1[np.argsort(ions1)[len(ions1)-top_n]] if len(ions1)>top_n else i>0
        ions0 = ions0[cutoff1]
        ions1 = ions1[cutoff1]
        ions = np.array([ions0,ions1])
    # Remove peaks outside the min_frag_mz to max_frag_mz range
    if min_frag_mz > 0:
        cutoff0 = ions0/max(ions0) >= min_frag_mz
        ions0 = ions0[cutoff0]
        ions1 = ions1[cutoff0]
    if max_frag_mz > 0:
        cutoff0 = ions0/max(ions0) <= max_frag_mz
        ions0 = ions0[cutoff0]
        ions1 = ions1[cutoff0]
    # # Duplicate m/z measurement
    check = len(np.unique(ions0)) != len(ions0)
    if check == True:
        temp = ions.copy()
        temp = pd.DataFrame(temp).T
        temp = temp[temp.groupby(0)[1].rank(ascending=False)<2]
        temp.drop_duplicates(subset=0, inplace=True)
        ions = np.array(temp.T)
    # Deisotope (experimental)
    if deiso: # TODO: Intensity, 2C13
        ions = deisotope(ions, m_proton, max_charge)
    return(ions)

def spscore(sub_spec, matched_ions, ftol, seq, mfrags):
    if mfrags.size > 0:
        # Consecutive fragments
        f = np.unique(mfrags)
        # B series
        bf = f[np.char.startswith(f, 'b')]
        if bf.size > 0:
            bf = np.array([i.replace('b' , '') for i in bf])
            bf3 = bf[np.char.endswith(bf, '+++')]
            bf = bf[np.invert(np.char.endswith(bf, '+++'))]
            bf2 = bf[np.char.endswith(bf, '++')]
            bf = bf[np.invert(np.char.endswith(bf, '++'))]
            bf3 = np.array([int(i.replace('+' , '')) for i in bf3])
            bf2 = np.array([int(i.replace('+' , '')) for i in bf2])
            bf = np.array([int(i.replace('+' , '')) for i in bf])
        else: 
            bf = bf2 = bf3 = np.array([0])
        # Y series
        yf = f[np.char.startswith(f, 'y')]
        if yf.size > 0:
            yf = np.array([i.replace('y' , '') for i in yf])
            yf3 = yf[np.char.endswith(yf, '+++')]
            yf = yf[np.invert(np.char.endswith(yf, '+++'))]
            yf2 = yf[np.char.endswith(yf, '++')]
            yf = yf[np.invert(np.char.endswith(yf, '++'))]
            yf3 = np.array([int(i.replace('+' , '')) for i in yf3])
            yf2 = np.array([int(i.replace('+' , '')) for i in yf2])
            yf = np.array([int(i.replace('+' , '')) for i in yf])
        else: 
            yf = yf2 = yf3 = np.array([0])
        # Claculate continuity per charge
        beta = (((sum((bf[1:]-bf[:-1])==1) + sum((yf[1:]-yf[:-1])==1)) * 0.075) +
                ((sum((bf2[1:]-bf2[:-1])==1) + sum((yf2[1:]-yf2[:-1])==1)) * 0.075) +
                ((sum((bf3[1:]-bf3[:-1])==1) + sum((yf3[1:]-yf3[:-1])==1)) * 0.075))
        # Immonium_ions
        rho = 0
        immonium_ions = [('H', 110.0718), ('Y', 136.0762), ('W', 159.0922), ('M', 104.0534), ('F', 120.0813)]
        immonium_ions = [(i[0],i[1]) for i in immonium_ions if i[1] >= sub_spec[0].min()-(i[1]/((1/ftol)*1E6))]
        for i in immonium_ions:
            if min(abs(sub_spec[0] - i[1])) <= (i[1]/((1/ftol)*1E6)): # immonium ion found
                minloc = np.where(abs(sub_spec[0]-i[1]) == min(abs(sub_spec[0] - i[1])))
                if sub_spec[1][minloc] > 0: # TODO minimum intensity to be considered?
                    if i[0] in seq: # increase rho if immonium in sequence
                        rho += 0.15
                    else: # decrease rho if immonium not in sequence
                        rho -= 0.15
        nm = len(mfrags)
        im = matched_ions
        nt = len(seq) * 2 * 3 # 2 series, 3 charge states # TODO param
        sp = (im * nm * (1 + beta) * (1 + rho)) / nt
    else:
        sp = 0
    return(sp)

def insertMods(peptide, mods):
    mods = mods.split(sep=", ")
    modlist = []
    for m in mods:
        #value = float(re.findall("\d+\.\d+", m)[0])
        pos, value = re.findall('[^()]+', m)
        value = float(value)
        if len(re.findall(r'\d+', pos)) > 0:
            pos = int(re.findall(r'\d+', pos)[0])
        elif pos == "N-term":
            pos = 1
        elif pos == "C-term":
            pos = len(peptide)
        modlist.append([pos, value])
    modlist = modlist[::-1] # Reverse list so it can be added without breaking pos
    omod = []
    opos = []
    for pos, value in modlist:
        peptide = peptide[:pos] + '[' + str(value) + ']' + peptide[pos:]
        omod.append(value)
        opos.append(pos-1)
    return(peptide, omod, opos)

def getTheoMH(sequence, nt, ct, mass,
              m_proton, m_hydrogen, m_oxygen):
    '''    
    Calculate theoretical MH using the PSM sequence.
    '''
    AAs = dict(mass._sections['Aminoacids'])
    MODs = dict(mass._sections['Fixed Modifications'])
    # total_aas = 2*m_hydrogen + m_oxygen
    total_aas = m_proton
    # total_aas += charge*m_proton
    #total_aas += float(MODs['nt']) + float(MODs['ct'])
    if nt:
        total_aas += float(MODs['nt'])
    if ct:
        total_aas += float(MODs['ct'])
    for i, aa in enumerate(sequence):
        if aa.lower() in AAs:
            total_aas += float(AAs[aa.lower()])
        if aa.lower() in MODs:
            total_aas += float(MODs[aa.lower()])
        # if aa.islower():
        #     total_aas += float(MODs['isolab'])
        # if i in pos:
        #     total_aas += float(mods[pos.index(i)]) TODO: add mod mass outside
    # MH = total_aas - m_proton
    return(total_aas)

def theoSpectrum(seq, blist, ylist, mods, pos, mass,
                 m_proton, m_hydrogen, m_oxygen, charge, dm=0):
    ## Y SERIES ##
    outy = []
    for i in ylist:
        yn = list(seq[-i:])
        if i < len(seq): nt = False
        else: nt = True
        fragy = getTheoMH(yn,nt,True,mass,
                          m_proton,m_hydrogen,m_oxygen) + 2*m_hydrogen + m_oxygen + dm
        outy += [fragy]
    ## B SERIES ##
    outb = []
    for i in blist:
        bn = list(seq[:i][::-1])
        if i > 0: ct = False
        else: ct = True
        fragb = getTheoMH(bn,True,ct,mass,
                          m_proton,m_hydrogen,m_oxygen) + dm # TODO only add +dm to fragments up until n_pos
        outb += [fragb]
    ## ADD FIXED MODS ## # TODO two modes, use mods from config file or input table
    # for i, m in enumerate(mods):
        # bpos = range(0, pos[mods.index(i)]+1)
        # ypos = range(len(seq)-pos[mods.index(i)]-1, len(seq))
        # bpos = pos[i]
        # ypos = len(seq)-pos[i]-1
        # spec[0] = spec[0][:bpos] + [b + m for b in spec[0][bpos:]]
        # spec[1] = spec[1][:ypos] + [y + m for y in spec[1][ypos:]]
    ## FRAGMENT MATRIX ##
    spec = []
    for c in range(1, charge+1):
        if c > 1:
            coutb = [(i+(c-1)*m_proton)/c for i in outb]
            couty = [(i+(c-1)*m_proton)/c for i in outy]
        else:
            coutb, couty = outb, outy
        spec += [[coutb, couty]]
    return(spec)

def addMod(spec, dm, pos, len_seq, blist, ylist):
    ## ADD MOD TO SITES ##
    bpos = [i >= pos+1 for i in blist]
    ypos = [i >= len_seq-pos for i in ylist][::-1]
    spec[0] = [spec[0][i]+dm if bpos[i]==True else spec[0][i] for i in list(range(0,len(spec[0])))]
    spec[1] = [spec[1][i]+dm if ypos[i]==True else spec[1][i] for i in list(range(0,len(spec[1])))]
    return(spec)

def makeFrags(seq, ch, full_y): # TODO: SLOW
    '''
    Name all fragments.
    '''
    bp = bh = yp = yh = []
    # if 'P' in seq[:-1]: # Cannot cut after
    #     # Disallowed b ions
    #     bp = [pos+1 for pos, char in enumerate(seq) if char == 'P']
    #     # Disallowed y ions
    #     yp = [pos for pos, char in enumerate(seq[::-1]) if char == 'P']
    # if 'H' in seq[1:]: # Cannot cut before
    #     # Disallowed b ions
    #     bh = [pos for pos, char in enumerate(seq) if char == 'H']
    #     # Disallowed y ions
    #     yh = [pos+1 for pos, char in enumerate(seq[::-1]) if char == 'H']
    seq_len = len(seq)
    blist = list(range(1,seq_len))
    blist = [i for i in blist if i not in bp + bh]
    ylist = list(range(1+int(full_y),seq_len+1))
    ylist = [i for i in ylist if i not in yp + yh]
    frags = []
    frags_m = []
    for c in range(1, ch+1):
        frags += [[["b" + str(i) + "+"*c for i in blist], ["y" + str(i) + "+"*c for i in ylist]]]
        frags_m += [[["b" + str(i) + "*" + "+"*c for i in blist], ["y" + str(i) + "*" + "+"*c for i in ylist]]]
    return(frags, frags_m, blist, ylist)

def fragCheck(plainseq, blist, ylist, dm_pos, charge):
    # ballowed = ['b'+str(i)+'*' if i >= dm_pos+1 else 'b'+str(i) for i in blist] * charge
    # yallowed = ['y'+str(i)+'*' if i >= len(plainseq)-dm_pos else 'y'+str(i) for i in ylist] * charge
    # cballowed = list(itertools.chain.from_iterable([['+'*i]*len(blist) for i in range(1,charge+1)]))
    # cyallowed = [['+'*i]*len(ylist) for i in range(1,charge+1)]
    if charge == 1:
        allowed = (['b'+str(i)+'*+' if i >= dm_pos+1 else 'b'+str(i)+'+' for i in blist] +
                   ['y'+str(i)+'*+' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+' for i in ylist])
    elif charge == 2:
        allowed = (['b'+str(i)+'*+' if i >= dm_pos+1 else 'b'+str(i)+'+' for i in blist] +
                   ['y'+str(i)+'*+' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+' for i in ylist] +
                   ['b'+str(i)+'*++' if i >= dm_pos+1 else 'b'+str(i)+'++' for i in blist] +
                   ['y'+str(i)+'*++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'++' for i in ylist])
    elif charge == 3:
        allowed = (['b'+str(i)+'*+' if i >= dm_pos+1 else 'b'+str(i)+'+' for i in blist] +
                   ['y'+str(i)+'*+' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+' for i in ylist] +
                   ['b'+str(i)+'*++' if i >= dm_pos+1 else 'b'+str(i)+'++' for i in blist] +
                   ['y'+str(i)+'*++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'++' for i in ylist] +
                   ['b'+str(i)+'*+++' if i >= dm_pos+1 else 'b'+str(i)+'+++' for i in blist] +
                   ['y'+str(i)+'*+++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+++' for i in ylist])
    else:
        allowed = (['b'+str(i)+'*+' if i >= dm_pos+1 else 'b'+str(i)+'+' for i in blist] +
                   ['y'+str(i)+'*+' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+' for i in ylist] +
                   ['b'+str(i)+'*++' if i >= dm_pos+1 else 'b'+str(i)+'++' for i in blist] +
                   ['y'+str(i)+'*++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'++' for i in ylist] +
                   ['b'+str(i)+'*+++' if i >= dm_pos+1 else 'b'+str(i)+'+++' for i in blist] +
                   ['y'+str(i)+'*+++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'+++' for i in ylist] +
                   ['b'+str(i)+'*++++' if i >= dm_pos+1 else 'b'+str(i)+'++++' for i in blist] +
                   ['y'+str(i)+'*++++' if i >= len(plainseq)-dm_pos else 'y'+str(i)+'++++' for i in ylist])
    return(allowed)

def findClosest(dm, dmdf, dmtol, pos):
    cand = [i for i in range(len(dmdf[1])) if dmdf[1][i] > dm-dmtol and dmdf[1][i] < dm+dmtol]
    closest = pd.DataFrame([dmdf[0][cand], dmdf[1][cand], dmdf[2][cand]]).T
    closest = closest.iloc[:, :3]
    closest.columns = ['name', 'mass', 'site']
    closest = pd.concat([closest, pd.Series({'name':'EXPERIMENTAL', 'mass':dm, 'site':[pos]}).to_frame().T], ignore_index=True)
    return(closest)

def findPos(dm_set, plainseq): # TODO fix sites now that this is array instead of DF
    def _where(sites, plainseq):
        sites = sites.site
        subpos = []
        for s in sites:
            if s == 'Anywhere' or s == 'exp':
                subpos = list(range(0, len(plainseq)))
                break
            elif s == 'Non-modified':
                subpos = [-1]
                break
            elif s == 'N-term':
                subpos = [0]
            elif s == 'C-term':
                subpos = [len(plainseq) - 1]
            else:
                subpos = subpos + list(np.where(np.array(list(plainseq)) == str(s))[0])
        subpos = list(dict.fromkeys(subpos))
        subpos.sort()
        return(subpos)
    dm_set['idx'] = dm_set.apply(_where, plainseq=plainseq, axis=1)
    dm_set = dm_set[dm_set.idx.apply(lambda x: len(x)) > 0]
    return(dm_set)

def getClosestIon(spectrum, ion):
    pos = bisect_left(spectrum, ion)
    if pos == 0:
        return spectrum[0]
    if pos == len(spectrum):
        return spectrum[-1]
    before = spectrum[pos - 1]
    after = spectrum[pos]
    if after - ion < ion - before:
        return(after)
    else:
        return(before)
    
def hyperscore(exp_spec, theo_spec, frags, ftol):
    assigned_mz = [getClosestIon(exp_spec[0], mz) for mz in theo_spec]
    assigned_ppm = np.absolute(np.divide(np.subtract(assigned_mz, theo_spec), theo_spec)*1000000)
    assigned_mask = assigned_ppm <= ftol
    assigned_mz = list(itertools.compress(assigned_mz, assigned_mask))
    if len(assigned_mz) == 0: return([], [], [], 0, 0, 0, 0)
    else:
        assigned_frags = list(itertools.compress(frags, assigned_mask))
        assigned_int = [exp_spec[1][list(exp_spec[0]).index(mz)] for mz in assigned_mz]
        assigned_int_mask = [f[0]=='b' for f in assigned_frags]
        i_b = sum(list(itertools.compress(assigned_int, assigned_int_mask)))
        i_y = sum(list(itertools.compress(assigned_int, ~np.array(assigned_int_mask))))
        i_sum = i_b + i_y
        n_b = len(set([f.replace('+', '') for f in assigned_frags if f[0]=='b']))
        n_y = len(set([f.replace('+', '') for f in assigned_frags if f[0]=='y']))
        if i_b == 0: i_b = 1
        if i_y == 0: i_y = 1
        hs = math.log((i_b) * (i_y)) + math.log(math.factorial((n_b))) + math.log(math.factorial(n_y))
    return(assigned_mz, assigned_int, assigned_frags, n_b, n_y, i_sum, hs)

def miniVseq(sub, plainseq, mods, pos, mass, ftol, dmtol, dmdf, m_proton, m_hydrogen, m_oxygen, score_mode, full_y):
    ## ASSIGNDB ##
    # assigndblist = []
    # assigndb = []
    ## FRAGMENT NAMES ##
    charge = sub.Charge
    if charge >= 4: charge = 4
    frags, frags_m, blist, ylist = makeFrags(plainseq, charge, full_y)
    ## DM ##
    exp_pos = 'exp'
    dm_set = findClosest(sub.DM, dmdf, dmtol, exp_pos) # Contains experimental DM
    dm_set = findPos(dm_set, plainseq)
    dm_set = dm_set[dm_set.mass != 0]
    theo_spec = theoSpectrum(plainseq, blist, ylist, mods, pos, mass,
                             m_proton, m_hydrogen, m_oxygen, charge)
    flat_theo_spec = sum(sum(theo_spec, []), [])
    flat_frags = sum(sum(frags, []), [])
    f_len = len(flat_frags)
    
    ## NON-MODIFIED ##
    NM_mz, NM_int, NM_frags, NM_n_b, NM_n_y, NM_i, NM_hs = hyperscore(sub.Spectrum, flat_theo_spec, flat_frags, ftol)
    nm_results = [NM_n_b+NM_n_y, NM_i, NM_hs, 'Non-modified', 0, None, NM_frags]
    
    ## DM OPERATIONS ##
    mod_results = [0, 0, 0, None, None, None, None]
    mod_site_range = []
    best_mod = 0
    hyb_results = [0, 0, 0, None, None, None, None]
    hyb_site_range = []
    best_hyb = 0
    exp_results = [0, 0, 0, None, None, None, None]
    exp_site_range = []
    place = 0
    for index, row in dm_set.iterrows():
        temp_hyb_site_range = []
        temp_mod_site_range = []
        dm = row.mass
        # TODO support both HYBRID and MOD scoring
        for dm_pos in row.idx:
            ## MOD HYPERSCORE ##
            allowed_mod = fragCheck(plainseq, blist, ylist, dm_pos, charge) # TODO support charge states > 4
            theo_spec_mod = [flat_theo_spec[i]+dm if '*' in allowed_mod[i] else flat_theo_spec[i] for i in range(0, f_len)]
            MOD_mz, MOD_int, MOD_frags, MOD_n_b, MOD_n_y, MOD_i, MOD_hs = hyperscore(sub.Spectrum, theo_spec_mod, allowed_mod, ftol)
            ## HYBRID HYPERSCORE ##
            HYB_frags = [i for i in MOD_frags if i not in NM_frags]
            if len(HYB_frags) == 0:
                HYB_int, HYB_frags, HYB_n_b, HYB_n_y, HYB_i, HYB_hs = NM_int, NM_frags, NM_n_b, NM_n_y, NM_i, NM_hs
            else:
                HYB_int = [i for i in MOD_int if i not in NM_int] + NM_int
                HYB_frags += NM_frags
                HYB_n_b = len(set([f.replace('+', '').replace('*', '') for f in HYB_frags if f[0]=='b']))
                HYB_n_y = len(set([f.replace('+', '').replace('*', '') for f in HYB_frags if f[0]=='y']))
                HYB_int_mask = [f[0]=='b' for f in HYB_frags]
                HYB_i_b = sum(list(itertools.compress(HYB_int, HYB_int_mask)))
                HYB_i_y = sum(list(itertools.compress(HYB_int, ~np.array(HYB_int_mask))))
                HYB_i = HYB_i_b + HYB_i_y
                if HYB_i_b == 0: HYB_i_b = 1
                if HYB_i_y == 0: HYB_i_y = 1
                HYB_hs = math.log((HYB_i_b) * (HYB_i_y)) + math.log(math.factorial((HYB_n_b))) + math.log(math.factorial(HYB_n_y))
            ## STORE RESULTS ##
            dm_pos = plainseq[dm_pos] + str(dm_pos+1)
            if row['name'] == 'EXPERIMENTAL':
                if score_mode:
                    exp_site_range += [HYB_hs]
                    if HYB_hs > exp_results[2]: # EXPERIMENTAL
                        exp_results = [HYB_n_b+HYB_n_y, HYB_i, HYB_hs, row['name'], dm, dm_pos, HYB_frags]
                if not score_mode:
                    exp_site_range += [MOD_hs]
                    if MOD_hs > exp_results[2]: # MOD
                        exp_results = [MOD_n_b+MOD_n_y, MOD_i, MOD_hs, row['name'], dm, dm_pos, MOD_frags]
            else:
                temp_hyb_site_range += [HYB_hs]
                temp_mod_site_range += [MOD_hs]
                if HYB_hs > hyb_results[2]: # TODO handle score ties
                    hyb_results = [HYB_n_b+HYB_n_y, HYB_i, HYB_hs, row['name'], dm, dm_pos, HYB_frags]
                    best_hyb = place
                if MOD_hs > mod_results[2]:
                    mod_results = [MOD_n_b+MOD_n_y, MOD_i, MOD_hs, row['name'], dm, dm_pos, MOD_frags]
                    best_mod = place
        place += 1
        mod_site_range += [temp_mod_site_range]
        hyb_site_range += [temp_hyb_site_range]
    # Pick site in the middle of the range
    exp_dm_pos = [i for i in range(0,len(exp_site_range)) if exp_site_range[i]==max(exp_site_range)]
    if len(exp_dm_pos)>0:
        exp_results[5] = plainseq[exp_dm_pos[len(exp_dm_pos)//2]]+str(exp_dm_pos[len(exp_dm_pos)//2]+1)
    mod_site_range = mod_site_range[best_mod]
    mod_dm_pos = [i for i in range(0,len(mod_site_range)) if mod_site_range[i]==max(mod_site_range)]
    if len(mod_dm_pos)>0:
            mod_results[5] = plainseq[mod_dm_pos[len(mod_dm_pos)//2]]+str(mod_dm_pos[len(mod_dm_pos)//2]+1)
    hyb_site_range = hyb_site_range[best_hyb]
    hyb_dm_pos = [i for i in range(0,len(hyb_site_range)) if hyb_site_range[i]==max(hyb_site_range)]
    if len(hyb_dm_pos)>0:
        hyb_results[5] = plainseq[hyb_dm_pos[len(hyb_dm_pos)//2]]+str(hyb_dm_pos[len(hyb_dm_pos)//2]+1)
    
    return(nm_results+[[]], exp_results+[exp_site_range], mod_results+[mod_site_range], hyb_results+[hyb_site_range])

def parallelFragging(query, parlist):
    m_proton = parlist[4]
    scan = query.scannum
    charge = query.charge
    MH = query.precursor_neutral_mass + (m_proton)
    MZ = (MH+m_proton*(charge-1))/charge
    plain_peptide = query.peptide
    if pd.isnull(query.modification_info):
        sequence = plain_peptide
        mod = []
        pos = []
    else:
        sequence, mod, pos = insertMods(plain_peptide, query.modification_info)
    if parlist[7]=="mzml":
        spectrum = parlist[9][np.where(np.array(parlist[8])==scan)[0][0]]
        score_mode = parlist[10]
        full_y = parlist[11]
        preference = parlist[12]
    else:
        spectrum = query.spectrum
        score_mode = parlist[8]
        full_y = parlist[9]
        preference = parlist[10]
    dm = query.massdiff
    # TODO use calc neutral mass?
    # Make a Vseq-style query
    sub = pd.Series([scan, charge, MH, sequence, spectrum, dm],
                    index = ["FirstScan", "Charge", "MH", "Sequence", "Spectrum", "DM"])
    nm_r, exp_r, mod_r, hyb_r = miniVseq(sub, plain_peptide, mod, pos,
                                         parlist[0], parlist[1], parlist[2],
                                         parlist[3], parlist[4], parlist[5], parlist[6],
                                         score_mode, full_y)
    if preference: # Preference NM > MOD > EXP
        if score_mode: # HYBRID
            best_r = [nm_r, hyb_r, exp_r][np.argmax([nm_r[2], hyb_r[2], exp_r[2]])]
        else: # MOD
            best_r = [nm_r, mod_r, exp_r][np.argmax([nm_r[2], mod_r[2], exp_r[2]])]
    else: # Preference NM > EXP > MOD
        if score_mode: # HYBRID
            best_r = [nm_r, exp_r, hyb_r][np.argmax([nm_r[2], exp_r[2], hyb_r[2]])]
        else: # MOD
            best_r = [nm_r, exp_r, mod_r][np.argmax([nm_r[2], exp_r[2], mod_r[2]])]
    try:
        if len(best_r[6]) == 0: sp = 0
        else:
            spfrags = np.array([i.replace('*', '') for i in best_r[6]])
            sp = spscore(sub.Spectrum, best_r[0], parlist[1], query.peptide, spfrags)
    except TypeError:
        best_r = [0, 0, 0, 0, 0, 0, 0, 0]
        sp = 0
    return([MH, MZ, dm, exp_r[0], exp_r[2], nm_r[0], nm_r[2], mod_r, hyb_r,
            best_r[7], best_r[4], best_r[5], sequence, best_r[0], best_r[1],
            best_r[2], best_r[3], sp])

def makeSummary(df, outpath, infile, raw, dmlist, startt, endt, decoy, protein):
    
    smods = df.REFRAG_name.value_counts()
    smods = smods[smods.index!='EXPERIMENTAL']
    lsmods = []
    for i in range(0, len(smods)):
        lsmods += [str(list(smods)[i]), '\t', str(list(smods.index)[i]), '\n']
    lsmods = ''.join(lsmods)
        
    summary = (
    "DATE\t{date}\n"
    "FILE\t{infile}\n"
    "RAW\t{raw}\n"
    "DMLIST\t{dmlist}\n"
    "SEARCH TIME\t{time}\n"
    "TOTAL PSMs\t{total}\n"
    "TARGET PSMs\t{target}\n"
    "REFRAGGED PSMs\t{refrag}\t({perc}% of total)\n"
    "REFRAGGED TARGET PSMs\t{refragt}\t({perct}% of targets)\n\n"
    "THEORETICAL MODIFICATIONS FREQUENCY\n\n"
    "{smods}").format(date=datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
    infile=str(infile),
    raw=str(raw),
    dmlist=str(dmlist),
    time=str(endt-startt),
    total=str(len(df)),
    target=str(len(df[~df[protein].str.startswith(decoy)])),
    refrag=str(len(df[df.REFRAG_name!='EXPERIMENTAL'])),
    perc=str(round(len(df[df.REFRAG_name!='EXPERIMENTAL'])/len(df)*100,2)),
    refragt=str(len(df[(df.REFRAG_name!='EXPERIMENTAL')&(~df[protein].str.startswith(decoy))])),
    perct=str(round(len(df[(df.REFRAG_name!='EXPERIMENTAL')&(~df[protein].str.startswith(decoy))])/len(df[~df[protein].str.startswith(decoy)])*100,2)),
    smods=lsmods)
    
    with open(outpath, 'w') as f:
        f.write(summary)
    
    return

def formatSiteRange(score_range, peptide):
    if not len(peptide) == len(score_range):
        return(peptide)
    max_score = max(score_range)
    site_range = ''
    for i in range(len(peptide)):
        if score_range[i] == max_score:
            site_range += peptide[i].lower()
        else:
            site_range += peptide[i]
    return(site_range)

def globalFDR(df, score_column, prot_column, decoy_prefix, filter_target=False, filter_fdr=0):
    fl = len(df)
    df['row_index'] = range(1, len(df) + 1)
    df['REFRAG_Label'] = df.apply(lambda x: 'Decoy' if (x[prot_column][0:len(decoy_prefix)]==decoy_prefix) else 'Target', axis = 1)
    df.sort_values(by=[score_column, 'REFRAG_Label'], inplace=True, ascending=False)
    df['Rank'] = df.groupby('REFRAG_Label').cumcount()+1
    df['Rank_T'] = np.where(df['REFRAG_Label']=='Target', df['Rank'], 0)
    df['Rank_T'] = df['Rank_T'].replace(0, np.nan).ffill()
    df['Rank_D'] = np.where(df['REFRAG_Label'] == 'Decoy', df['Rank'], 0)
    df['Rank_D'] =  df['Rank_D'].replace(0, np.nan).ffill()
    df['REFRAG_FDR'] = df['Rank_D']/df['Rank_T']
    df['REFRAG_FDR'] =  df['REFRAG_FDR'].replace(np.nan, 0).ffill()
    df.sort_values(by='row_index', inplace=True, ascending=True)
    if filter_fdr > 0:
        logging.info("Filtering at " + str(filter_fdr*100) + "% FDR...")
        df = df[df.REFRAG_FDR <= filter_fdr]
    if filter_target:
        logging.info("Removing decoys...")
        df = df[~df.REFRAG_Label.str.startswith(decoy_prefix)]
    if filter_fdr > 0 or filter_target:
        logging.info(str(len(df)) + "out of " + str(fl) + " PSMs passed the filter(s).")
    df.drop(['row_index', 'Rank', 'Rank_T', 'Rank_D'], axis = 1, inplace = True)
    return(df)

def main(args):
    '''
    Main function
    '''
    # Parameters
    chunks = int(mass._sections['Search']['batch_size'])
    ftol = float(mass._sections['Search']['f_tol'])
    dmtol = float(mass._sections['Search']['dm_tol'])
    score_mode = bool(int(mass._sections['Search']['score_mode']))
    full_y = bool(int(mass._sections['Search']['full_y']))
    preference = bool(int(mass._sections['Search']['preference']))
    top_n = int(mass._sections['Spectrum Processing']['top_n'])
    bin_top_n = bool(int(mass._sections['Spectrum Processing']['bin_top_n']))
    min_ratio = float(mass._sections['Spectrum Processing']['min_ratio'])
    min_frag_mz = float(mass._sections['Spectrum Processing']['min_fragment_mz'])
    max_frag_mz = float(mass._sections['Spectrum Processing']['max_fragment_mz'])
    deiso = bool(int(mass._sections['Spectrum Processing']['deisotope']))
    decoy_prefix = str(mass._sections['FDR']['decoy_prefix'])
    prot_column = str(mass._sections['FDR']['prot_column'])
    filter_target =  bool(int(mass._sections['FDR']['filter_target']))
    filter_fdr = float(mass._sections['FDR']['filter_fdr'])
    m_proton = mass.getfloat('Masses', 'm_proton')
    m_hydrogen = mass.getfloat('Masses', 'm_hydrogen')
    m_oxygen = mass.getfloat('Masses', 'm_oxygen')
    # Debug
    debug_scores = bool(int(mass._sections['Debug']['debug_scores']))
    
    infiles, rawfiles, rawbase = preProcess(args)
    
    checked = checkParams(mass, infiles)
    if checked != 0:
        sys.exit("ERROR: Invalid parameters.")
        
    for f in infiles:
        if os.path.basename(f).split(sep=".")[0] in rawbase:
            infile = f
            rawfile = rawfiles[rawbase.index(os.path.basename(f).split(sep=".")[0])]
        else:
            logging.info("ERROR: No MS data file (mzML or MGF) found for MSFragger file " + str(f) + "\nSkipping...")
            continue

        # Read results file from MSFragger
        if len(args.dia) > 0:
            logging.info("Reading MSFragger file (" + str(os.path.basename(Path(infile[0:-8] + infile[-4:]))) + ")...")
            try: df = pd.read_csv(Path(infile[0:-8] + infile[-4:]), sep="\t")
            except pd.errors.ParserError:
                sep = '\t'
                coln = pd.read_csv(Path(infile[0:-8] + infile[-4:]), sep=sep, nrows=1).columns.tolist()
                max_columns = max(open(Path(infile[0:-8] + infile[-4:]), 'r'), key = lambda x: x.count(sep)).count(sep) + 1
                coln += ['extra_column_' + str(i) for i in range(0,(max_columns-len(coln)))]
                df = pd.read_csv(Path(infile[0:-8] + infile[-4:]), header=None, sep=sep, skiprows=1, names=coln)
        else:
            logging.info("Reading MSFragger file (" + str(os.path.basename(Path(infile))) + ")...")
            try: df = pd.read_csv(Path(infile), sep="\t")
            except pd.errors.ParserError:
                sep = '\t'
                coln = pd.read_csv(Path(infile), sep=sep, nrows=1).columns.tolist()
                max_columns = max(open(Path(infile), 'r'), key = lambda x: x.count(sep)).count(sep) + 1
                coln += ['extra_column_' + str(i) for i in range(0,(max_columns-len(coln)))]
                df = pd.read_csv(Path(infile), header=None, sep=sep, skiprows=1, names=coln)
        if len(args.scanrange) > 0:
            logging.info("\tFiltering scan range " + str(args.scanrange[0]) + "-" + str(args.scanrange[1]) + "...")
            df = df[(df.scannum>=args.scanrange[0])&(df.scannum<=args.scanrange[1])]
        logging.info("\t" + str(len(df)) + " lines read.")
        if len(df) < 1:
            logging.error("No scans to search. Stopping.")
            sys.exit()
        # Read raw file
        msdata, mode, index2 = readRaw(Path(rawfile))
        # Read DM file
        logging.info("Reading DM file (" + str(os.path.basename(Path(args.dmfile))) + ")...")
        dmdf = pd.read_csv(Path(args.dmfile), sep="\t") # TODO check for duplicates (when both DM and SITE is the same)
        dmdf = dmdf.iloc[:, :3]
        dmdf.columns = ["name", "mass", "site"]
        dmdf.site = dmdf.site.apply(literal_eval)
        dmdf.site = dmdf.apply(lambda x: list(dict.fromkeys(x.site)), axis=1)
        dmdf = dmdf.T.to_numpy()
        #dmdf[3] = np.array([''.join(set(''.join(literal_eval(i)).replace('N-term', '0').replace('C-term', '1'))) for i in dmdf[3]]) # Nt = 0, Ct = 1
        logging.info("\t" + str(len(dmdf[0])) + " theoretical DMs read.")
        # Prepare to parallelize
        logging.info("Refragging...")
        logging.info("\t" + "Locating scans...")
        starttime = datetime.now()
        if mode == "mzml":
            spectra = msdata.getSpectra()
            spectra_n = np.array([int(s.getNativeID().split("=")[-1]) for s in spectra])
            # if len(args.scanrange) > 0:
            #     logging.info("Filtering scan range " + str(args.scanrange[0]) + "-" + str(args.scanrange[1]) + "...")
            #     spectra = spectra[np.where((np.array(spectra_n)>=0)&(np.array(spectra_n)<=args.scanrange[0]))[0][0]:np.where((np.array(spectra_n)>=0)&(np.array(spectra_n)<=args.scanrange[-1]))[0][-1]]
            #     spectra_n = spectra_n[np.where((np.array(spectra_n)>=0)&(np.array(spectra_n)<=args.scanrange[0]))[0][0]:np.where((np.array(spectra_n)>=0)&(np.array(spectra_n)<=args.scanrange[-1]))[0][-1]]
            peaks = [s.get_peaks() for s in spectra]
            empty_peaks = set([i for i, p in enumerate(peaks) if len(p[0]) == 0])
            # Skip empty spectra
            if len(empty_peaks) > 0:
                peaks = [p for i, p in enumerate(peaks) if i not in empty_peaks]
                spectra_n = [s for i, s in enumerate(spectra_n) if i not in empty_peaks]
            ions = [np.array([p[0], p[1]]) for p in peaks]
            ions0 = [np.array(p[0]) for p in peaks]
            ions1 = [np.array(p[1]) for p in peaks]
            # Normalize intensity
            logging.info("\t" + "Normalise intensity...")
            ions1 = [(ions1[i]/max(ions1[i]))*100 for i in range(len(ions))]
            # Remove peaks below min_ratio
            logging.info("\t" + "Filter by ratio...")
            if min_ratio > 0:
                cutoff1 = [i/max(i) >= min_ratio for i in ions1]
                ions0 = [ions0[i][cutoff1[i]] for i in range(len(ions))]
                ions1 = [ions1[i][cutoff1[i]] for i in range(len(ions))]
            # Remove peaks outside the min_frag_mz to max_frag_mz range
            if min_frag_mz > 0:
                cutoff0 = ions0/max(ions0) >= min_frag_mz
                ions0 = ions0[cutoff0]
                ions1 = ions1[cutoff0]
            if max_frag_mz > 0:
                cutoff0 = ions0/max(ions0) <= max_frag_mz
                ions0 = ions0[cutoff0]
                ions1 = ions1[cutoff0]
            # Return only top N peaks
            logging.info("\t" + "Filter top N...")
            if top_n > 0:
                cutoff1 = [i >= i[np.argsort(i)[len(i)-top_n]] if len(i)>top_n else i>0 for i in ions1]
                ions0 = [ions0[i][cutoff1[i]] for i in range(len(ions))]
                ions1 = [ions1[i][cutoff1[i]] for i in range(len(ions))]
            ions = [np.array([ions0[i],ions1[i]]) for i in range(len(ions))]
            # # Duplicate m/z measurement
            logging.info("\t" + "Filter duplicate m/z measurements...")
            check = [len(np.unique(i)) != len(i) for i in ions0]
            for i in range(len(check)):
                if check[i] == True:
                    temp = ions[i].copy()
                    temp = pd.DataFrame(temp).T
                    temp = temp[temp.groupby(0)[1].rank(ascending=False)<2]
                    temp.drop_duplicates(subset=0, inplace=True)
                    ions[i] = np.array(temp.T)
            cc = 0
            for i in ions:
                if len(i[0])<2:
                    cc += 1
            logging.info("\t" + str(cc) + " empty spectra...")
        else:
            df.scannum = df.scannum.astype(int)
            df["spectrum"] = df.apply(lambda x: locateScan(x.scannum, mode, msdata,
                                                           0, 0, index2, top_n,
                                                           bin_top_n, min_ratio,
                                                           min_frag_mz, max_frag_mz,
                                                           m_proton, deiso, x.charge),
                                      axis=1)
        indices, rowSeries = zip(*df.iterrows())
        rowSeries = list(rowSeries)
        tqdm.pandas(position=0, leave=True)
        if len(df) <= chunks:
            chunks = math.ceil(len(df)/args.n_workers)
        if mode == "mzml":
            parlist = [mass, ftol, dmtol, dmdf, m_proton, m_hydrogen, m_oxygen, mode, spectra_n, ions, score_mode, full_y, preference]
        else:
            parlist = [mass, ftol, dmtol, dmdf, m_proton, m_hydrogen, m_oxygen, mode, score_mode, full_y, preference]
        logging.info("\tBatch size: " + str(chunks) + " (" + str(math.ceil(len(df)/chunks)) + " batches)")
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:
            refrags = list(tqdm(executor.map(parallelFragging,
                                             rowSeries,
                                             itertools.repeat(parlist),
                                             chunksize=chunks),
                                total=len(rowSeries)))
        if 'spectrum' in df.columns:
            df = df.drop('spectrum', axis = 1)
        df['templist'] = refrags
        # Experimental information #
        df['REFRAG_MH'] = pd.DataFrame(df.templist.tolist()).iloc[:, 0]. tolist()
        df['REFRAG_exp_MZ'] = pd.DataFrame(df.templist.tolist()).iloc[:, 1]. tolist()
        df['REFRAG_exp_DM'] = pd.DataFrame(df.templist.tolist()).iloc[:, 2]. tolist()
        df['REFRAG_exp_ions_matched'] = pd.DataFrame(df.templist.tolist()).iloc[:, 3]. tolist()
        df['REFRAG_exp_hyperscore'] = pd.DataFrame(df.templist.tolist()).iloc[:, 4]. tolist()
          # TODO: df['REFRAG_exp_hyperscore_M']
          # TODO: df['REFRAG_exp_hyperscore_H']
        # Non-modified information #
        df['REFRAG_nm_ions_matched'] = pd.DataFrame(df.templist.tolist()).iloc[:, 5]. tolist()
        df['REFRAG_nm_hyperscore'] = pd.DataFrame(df.templist.tolist()).iloc[:, 6]. tolist()
        # Modified information # (TESTING)
        if debug_scores:
            df['REFRAG_MOD_hs'] = pd.DataFrame(df.templist.tolist()).iloc[:, 7]. tolist()
            df['REFRAG_HYB_hs'] = pd.DataFrame(df.templist.tolist()).iloc[:, 8]. tolist()
        # Best candidate information #
        df['REFRAG_score_range'] = pd.DataFrame(df.templist.tolist()).iloc[:, 9]. tolist()
        df['REFRAG_site_range'] = df.apply(lambda x: formatSiteRange(x.REFRAG_score_range, x.peptide), axis=1)
        df['REFRAG_DM'] =  pd.DataFrame(df.templist.tolist()).iloc[:, 10]. tolist()
        df['REFRAG_site'] = pd.DataFrame(df.templist.tolist()).iloc[:, 11]. tolist()
        df['REFRAG_sequence'] = pd.DataFrame(df.templist.tolist()).iloc[:, 12]. tolist() # TODO: add DM?
        df['REFRAG_ions_matched'] = pd.DataFrame(df.templist.tolist()).iloc[:, 13]. tolist()
        df['REFRAG_sum_intensity'] = pd.DataFrame(df.templist.tolist()).iloc[:, 14]. tolist()
        df['REFRAG_hyperscore'] = pd.DataFrame(df.templist.tolist()).iloc[:, 15]. tolist()
        df['REFRAG_name'] = pd.DataFrame(df.templist.tolist()).iloc[:, 16]. tolist()
        df['REFRAG_sp_score'] = pd.DataFrame(df.templist.tolist()).iloc[:, 17]. tolist()
        df = df.drop('templist', axis = 1)
        logging.info("Calculating FDR...")
        df = globalFDR(df, 'REFRAG_hyperscore', prot_column, decoy_prefix, filter_target, filter_fdr)
        
        try:
            refragged = len(df)-df.REFRAG_name.value_counts()['EXPERIMENTAL']
        except KeyError:
            refragged = len(df)
        prefragged = round((refragged/len(df))*100,2)
        logging.info("\t" + str(refragged) + " (" + str(prefragged) + "%) refragged PSMs.")
        endtime = datetime.now()

        logging.info("Preparing workspace...")
        # get the name of script
        script_name = os.path.splitext( os.path.basename(__file__) )[0].upper()
        # if output directory is not defined, get the folder from given file + script name
        outdir = args.outdir if args.outdir else os.path.join(os.path.dirname(infile), script_name.lower())
        os.makedirs(outdir, exist_ok=True)
        basename = os.path.basename(infile)
        outpath = Path(outdir) / basename
        outsum = Path(outdir) / f"{Path(basename).stem}_SUMMARY{Path(basename).suffix}"

        if len(args.scanrange) > 0:
            outpath = Path(outdir) / f"{outpath.stem}_{str(args.scanrange[0]) + '-' + str(args.scanrange[1])}{outpath.suffix}"
        logging.info(f"Writing output file {outpath}...")
        df.to_csv(outpath, index=False, sep='\t', encoding='utf-8')
        logging.info("Done.")

        if len(args.scanrange) > 0:
            outsum = Path(outdir) / f"{outsum.stem}_{str(args.scanrange[0]) + '-' + str(args.scanrange[1])}{outsum.suffix}"
        logging.info(f"Writing summary file {outsum}...")
        makeSummary(df, outsum, infile, rawfile, args.dmfile, starttime, endtime, decoy_prefix, prot_column)
        logging.info("Done.")

    return

if __name__ == '__main__':

    # multiprocessing.freeze_support()
    # parse arguments
    parser = argparse.ArgumentParser(
        description='ReFrag',
        epilog='''
        Example:
            python ReFrag.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "ReFrag.ini")
    
    # TODO parameter: exclude DM range (consider as NM)= default (-3, 0)
    parser.add_argument('-i',  '--infile', required=True, help='MSFragger results file')
    parser.add_argument('-r',  '--rawfile', required=True, help='MS Data file (MGF or MZML)')
    parser.add_argument('-d',  '--dmfile', required=True, help='DeltaMass file')
    parser.add_argument('-a', '--dia', default=[], help='DIA mode (looks for MS Data files with _chN suffix)',
                        type=lambda s: [int(item) for item in s.split(',')])
    parser.add_argument('-s', '--scanrange', default=[], help='Scan range to search',
                        type=lambda s: [int(item) for item in s.split(',')])
    parser.add_argument('-o', '--outdir', help='Path to the output directory')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    parser.add_argument('-w',  '--n_workers', type=int, default=os.cpu_count(), help='Number of threads/n_workers')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    mass = configparser.ConfigParser(inline_comment_prefixes='#')
    if Path(args.config).is_file():
        mass.read(args.config)
    else:
        sys.exit('ERROR: Cannot find configuration file at' + str(args.config))
    # if something is changed, write a copy of ini
    if mass.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/RefMod.ini', 'w') as newconfig:
            mass.write(newconfig)

    # logging debug level. By default, info level
    if os.path.isdir(args.infile):
        log_file = os.path.join(args.infile + '/RefMod.log')
        log_file_debug = os.path.join(args.infile + '/RefMod_debug.log')
    elif '*' in args.infile:
        log_file = os.path.join(os.path.dirname(args.infile) + '/' + os.path.basename(args.infile).replace("*", "") + '_RefMod.log')
        log_file_debug = os.path.join(os.path.dirname(args.infile) + '/' + os.path.basename(args.infile).replace("*", "") + '_RefMod_debug.log')
    else:
        log_file = args.infile[:-4] + '_RefMod.log'
        log_file_debug = args.infile[:-4] + '_RefMod_debug.log'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])
        
    if len(args.scanrange) > 2:
        sys.exit('ERROR: Scan range list can only contain two elements (start and end)')

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')