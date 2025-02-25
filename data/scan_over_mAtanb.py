import sys
import csv
import subprocess
import pandas as pd
import pyslha
import numpy as np
import math
import random
from multiprocessing import Pool

csvfile = pd.read_csv(sys.argv[1], header=None, delim_whitespace=True)

input_file = """
Block MODSEL                 # Select model
   12    125                # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   4                    # pole mass loop order
    5   1                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   3                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L, 3 = 4L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   0                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   0                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   1                    # loop library (0 = softsusy)
   32   2                    # loop level to calculate AMM
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.20               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   172.57      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
     1     1.85000000E+00   # TanBeta
     2     1.00000000E+03   # mAIN
     3     4.71405000E-01   # ptIN
     4    -7.47525000E-01   # p1IN
     5    -1.16169000E+00   # p2IN
     6    -9.69514000E-01   # p3IN
     7     0.00000000E+00   # p4IN
     8     0.00000000E+00   # p5IN
     9     0.00000000E+00   # p6IN
    10     0.00000000E+00   # p7IN
    11    -1.19492000E-01   # qtIN
    12     2.73834000E+00   # q1IN
    13     3.77176000E+00   # q2IN
    14     3.34711000E+00   # q3IN
    15     0.00000000E+00   # q4IN
    16     0.00000000E+00   # q5IN
    17     0.00000000E+00   # q6IN
    18     0.00000000E+00   # q7IN
    19     1.22814000E+00   # rtIN
    20     5.85514000E+00   # r1IN
    21    -4.10998000E+00   # r2IN
    22     6.09427000E+00   # r3IN
    23    -7.86293000E+00   # r4IN
    24     0.00000000E+00   # r5IN
    25     0.00000000E+00   # r6IN
    26     0.00000000E+00   # r7IN
Block EXTPAR
    0   1e+7                  # Qin
    1   125.
Block FlexibleDecay
   0   1                    # calculate decays (0 = no, 1 = yes)
   1   1e-5                 # minimum BR to print
   2   4                    # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
   3   1                    # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
   4   2                    # off-shell decays into VV pair
   5   1
   6   1
   7   1
   8   0
   9   1
"""

def get_spectrum(x):
    slha_input = pyslha.readSLHA(input_file, ignorenomass=True)

    window = 0.05
    slha_input.blocks['MINPAR'][1] = random.uniform(1.8, 1.92)
    slha_input.blocks['MINPAR'][2] = random.uniform(500, 2500)

    slhastring = pyslha.writeSLHA(slha_input, ignorenobr=True, precision=16)
    proc = subprocess.run(
             [
                "/home/wojciech/Programowanie/c++/FlexibleDevelopment/models/THDMIIRoC/run_THDMIIRoC.x",
                "--higgsbounds-dataset=/home/wojciech/HEP-software/HBDataSet",
                "--higgssignals-dataset=/home/wojciech/HEP-software/HSDataSet",
                #"--lilith-db=$HOME/HEP-software/Lilith/data/latestRun2.list",
                '--slha-input-file=-'
             ],
             stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=slhastring.encode()
    )
    try:
       spectrum = pyslha.readSLHA(proc.stdout.decode())
    except Exception as e:
       return [None, None]

    print(proc.stdout.decode())

    hb_r = 0.
    for _idx, val in spectrum.blocks["HIGGSBOUNDS"].items():
       if _idx[1] == 1:
          hb_r = max(val, hb_r)
    y = []
    for i in range(20):
        y.append(slha_input.blocks['MINPAR'][i+1])
    return [spectrum.blocks["HIGGSSIGNALS"][4], hb_r] + [spectrum.blocks["MASS"][6]] + [spectrum.blocks["MASS"][25]] + [spectrum.blocks["MINPAR"][1]] + [spectrum.blocks["MINPAR"][2]] + y

def f():
   return get_spectrum(csvfile.to_numpy()[398])

data = np.array(Pool(int(sys.argv[2])).starmap(f, [() for _ in range(int(sys.argv[3]))]))

with open(f'scan_over_mAtanb.csv', 'a', encoding='UTF8', newline='') as file:
   writer = csv.writer(file)
   writer.writerows(data)
