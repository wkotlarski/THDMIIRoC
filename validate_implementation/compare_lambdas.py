import sys
import subprocess
import pandas as pd
import pyslha
import matplotlib.pyplot as plt
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
    4   0                    # pole mass loop order
    5   0                    # EWSB loop order
    6   1                    # beta-functions loop order
    7   0                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   1                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
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
   24   0                    # individual threshold correction loop orders
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
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733000000e+02      # mtop(pole)
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
     1     0.2246540E+01   # TanBeta
     2     800   # mA
     3     0.4714050E+00   # ptIN
     4     -0.1166670E+01   # p1IN
     5     -0.1377120E+01  # p2IN
     6     0.00000000E+00   # p3IN
     7     0.00000000E+00   # p4IN
     8     0.00000000E+00   # p5IN
     9     0.00000000E+00   # p6IN
    10     0.00000000E+00   # p7IN
    11     -0.1194920E+00   # qtIN
    12     0.4198090E+01   # q1IN
    13     0.4606740E+01    # q2IN
    14     -0.1156890E+00   # q3IN
    15     0   # q4IN
    16     0   # q5IN
    17     0.00000000E+00   # q6IN
    18     0.00000000E+00   # q7IN
    19     0.1228140E+01   # rtIN
    20     0.9390920E+01   # r1IN
    21     0.1027310E+02   # r2IN
    22     -0.9635690E-02   # r3IN
    23     0.    # r4IN
    24     0.00000000E+00   # r5IN
    25     0.00000000E+00   # r6IN
    26     0.00000000E+00   # r7IN
Block EXTPAR
    0   1e+7                  # Qin
    1   125.
Block FlexibleDecay
   0   0                    # calculate decays (0 = no, 1 = yes)
   1   1e-5                 # minimum BR to print
   2   4                    # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
   3   1                    # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
   4   2                    # off-shell decays into VV pair
   5   1
   6   1
   7   1
   8   1
   9   1
"""

worst_idx = -1
worst_diff = 0
hs = -1
best_idx = -1
best_pval = 0
def get_spectrum(x, idx):
    global worst_diff
    global worst_idx
    global hs
    global best_idx
    global best_pval

    #if any(x==idx for x in [266,267,275,276,284,285,293,294,302,303,311,317,318,326,327]):
    #    return

    slha_input = pyslha.readSLHA(input_file, ignorenomass=True)

    offset = 0 #3
    slha_input.blocks['MINPAR'][1] = x[4]
    slha_input.blocks['MINPAR'][3] = x[18-offset]
    for i in range(7):
       slha_input.blocks['MINPAR'][4+i] = x[11+i-offset]
    slha_input.blocks['MINPAR'][11] = x[26-offset]
    for i in range(7):
       slha_input.blocks['MINPAR'][12+i] = x[19+i-offset]
    slha_input.blocks['MINPAR'][19] = x[34-offset]
    for i in range(7):
       slha_input.blocks['MINPAR'][20+i] = x[27+i-offset]

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
    except:
       return [None, None]

    print(proc.stdout.decode())
    print(f'index: {idx}')
    print(f'mhG: {x[6-offset]}, mhW: {spectrum.blocks["MASS"][25]}, diff (%): {100*(1-x[51-offset]/spectrum.blocks["MASS"][25])}')

    print(f'l1 high GP: {x[43-offset]}')
    print(f'l2 high GP: {x[44-offset]}')
    print(f'l3 high GP: {x[45-offset]}')
    print(f'l4 high GP: {x[46-offset]}')
    print(f'l5 high GP: {x[47-offset]}')
    print(f'Yt high GP: {x[50-offset]}')
    print(f'l1 low GP: {x[35-offset]}')
    print(f'l2 low GP: {x[36-offset]}')
    print(f'l3 low GP: {x[37-offset]}')
    print(f'l4 low GP: {x[38-offset]}')
    print(f'l5 low GP: {x[39-offset]}')
    print(f'Yt low GP: {x[40-offset]}')
    print(f'l1: {x[35-offset]/spectrum.blocks["HMIX"][31]}')
    print(f'l2: {x[36-offset]/spectrum.blocks["HMIX"][32]}')
    print(f'l3: {x[37-offset]/spectrum.blocks["HMIX"][33]}')
    print(f'l4: {x[38-offset]/spectrum.blocks["HMIX"][34]}')
    print(f'l5: {x[39-offset]/spectrum.blocks["HMIX"][35]}')
    print(f'top mass: {x[6]}')
    print(f'Higgs mass: {x[6]}')
    print(f'Higgs mass 2: {x[51]}')
    print(f'boundary mass: {x[5]}')
    if max(abs(2-x[35-offset]/spectrum.blocks["HMIX"][31]), abs(2-x[36-offset]/spectrum.blocks["HMIX"][32]), abs(1-x[37-offset]/spectrum.blocks["HMIX"][33]), abs(1-x[38-offset]/spectrum.blocks["HMIX"][34]), abs(1-x[39-offset]/spectrum.blocks["HMIX"][35]) if not math.isnan(x[39-offset]/spectrum.blocks["HMIX"][35]) else 0.) > worst_diff:
        worst_diff = max(abs(2-x[35-offset]/spectrum.blocks["HMIX"][31]), abs(2-x[36-offset]/spectrum.blocks["HMIX"][32]), abs(1-x[37-offset]/spectrum.blocks["HMIX"][33]), abs(1-x[38-offset]/spectrum.blocks["HMIX"][34]), abs(1-x[39-offset]/spectrum.blocks["HMIX"][35]) if not math.isnan(x[39-offset]/spectrum.blocks["HMIX"][35]) else 0.)
        worst_idx = idx

for idx, row in enumerate(csvfile.to_numpy()):
   get_spectrum(row, idx)

print('worst point')
print(f'max diff: {worst_diff}')
print(f'ind {worst_idx}')
get_spectrum(csvfile.to_numpy()[worst_idx], 1)
