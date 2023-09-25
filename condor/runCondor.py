import time
import os
import sys
import re
from argparse import ArgumentParser
from datetime import datetime

# Ref: https://www.geeksforgeeks.org/how-to-search-and-replace-text-in-a-file-in-python/
def replacetext(filename, search_text, replace_text):
    with open(filename, 'r+') as f:
        file = f.read()
        file = re.sub(search_text, replace_text, file)
        f.seek(0)
        f.write(file)
        f.truncate()
    # return 'Text {} replaced by {}'.format(search_text, replace_text)


os.makedirs('./condor_out/', exist_ok=True)
os.makedirs('./condor_out/out/', exist_ok=True)
os.makedirs('./condor_out/err/', exist_ok=True)
os.makedirs('./condor_out/log/', exist_ok=True)
os.system('rm condor_*.job')
os.system('rm ./condor_out/out/*')
os.system('rm ./condor_out/err/*')
os.system('rm ./condor_out/log/*')

# Usage: ./TrackletAna [NevtToRun] [skip] [layer] [randhit_case] [clusplit_case] [misalignment_num] [dRCut]
# Example: ./TrackletAna 2000 0 12 0 0 0 0.5
nEvtsPerRun = 2000
layers = [12, 23, 13]
skips = [0]
randhitcases = [1, 2, 3]
clustercases = [1, 2, 3]
misliannums = [100, 101]
# for i in range(0, 40):
#     misliannums.append((i*10) + 1)

drcut = [0.4, 0.6, 999] # For simple pion, 999 means no dR cut
drcuttext = ['0p4','0p6', '999']
# drcut = [0.4, 0.6] # For nominal analysis
# drcuttext = ['0p4','0p6']

# nominal run arguments
run_args = []
# for l in layers:
#     for iskip in skips: # skip
#         # arg = [[NevtToRun], [skip], [layer], [randhit_case], [clusplit_case], [misalignment_num], [dRCut]]
#         arg = [nEvtsPerRun, iskip, l, 0, 0, 0, 0.5, '0p5']
#         run_args.append(arg)

# Systematic uncertainty: random hit variations
# for l in layers:
#     for iskip in skips: # skip
#         for randhitcase in randhitcases:
#             arg = [nEvtsPerRun, iskip, l, randhitcase, 0, 0, 0.5, '0p5']
#             run_args.append(arg)

# Sytematic uncertainty: cluster splitting variations
# for l in layers:
#     for iskip in skips: # skip
#         for clustercase in clustercases:
#             arg = [nEvtsPerRun, iskip, l, 0, clustercase, 0, 0.5, '0p5']
#             run_args.append(arg)

# Misalignment study 
for l in layers:
    for iskip in skips: # skip
        for misalignnum in misliannums:
            arg = [nEvtsPerRun, iskip, l, 0, 0, misalignnum, 0.5, '0p5']
            run_args.append(arg)

# Systematic uncertainty: dR cut variations
# for l in layers:
#     for iskip in skips: # skip
#         for idr in range(len(drcut)):
#             arg = [nEvtsPerRun, iskip, l, 0, 0, 0, drcut[idr], drcuttext[idr]]
#             run_args.append(arg)

print('Total number of runs: {}'.format(len(run_args)))

for iarg in range(len(run_args)):
    arg = run_args[iarg]
    print('Run set {}, events Per Run = {}, skipping = {} events, layer = {}, random hit case = {}, cluster split case = {}, Mis-alignment num = {}, dR cut for tracklet reco = {}'.format(iarg, arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[6]))
    newfile = 'condor_Nevt{}_Skip{}_layer{}_RandhitCase{}_ClusterSplitCase{}_MisAlignNum{}_dr{}.job'.format(arg[0], arg[1], arg[2], arg[3], arg[4], arg[5], arg[7])
    os.system('cp condor.job {}'.format(newfile))
    replacetext(newfile, 'NEVENTS', '{}'.format(arg[0]))
    replacetext(newfile, 'SKIP', '{}'.format(arg[1]))
    replacetext(newfile, 'LAYER', '{}'.format(arg[2]))
    replacetext(newfile, 'RANDHITCASE', '{}'.format(arg[3]))
    replacetext(newfile, 'CLUSTERSPLITCASE', '{}'.format(arg[4]))
    replacetext(newfile, 'MISALIGNNUM', '{}'.format(arg[5]))
    replacetext(newfile, 'DRCUT', '{}'.format(arg[6]))
    cmd = 'condor_submit ' + newfile
    os.system(cmd)