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

layers = [12, 23, 13]
# misliannums = [(i*10) + 1 for i in range(0, 40)]
misliannums = [0]
for i in range(0, 40):
    misliannums.append((i*10) + 1)

nEvtsPerRun = 4000
for l in layers:
    for i in range(0, 1):
        for randhitcase in range(0, 4):
            for misalignnum in misliannums: 
                if randhitcase > 0 and misalignnum > 0:
                    continue
                
                print('Run set {}, skipping {} events, layer {}, random hit case {}, Mis-alignment num {}'.format(i, i*nEvtsPerRun, l, randhitcase, misalignnum))
                newfile = 'condor_layer{}_Nevt{}_Skip{}_RandhitCase{}_MisAlignNum{}.job'.format(l, nEvtsPerRun, i*nEvtsPerRun, randhitcase, misalignnum)
                os.system('cp condor.job {}'.format(newfile))
                replacetext(newfile, 'NEVENTS', '{}'.format(nEvtsPerRun))
                replacetext(newfile, 'SKIP', '{}'.format(i*nEvtsPerRun))
                replacetext(newfile, 'LAYER', '{}'.format(l))
                replacetext(newfile, 'RANDHITCASE', '{}'.format(randhitcase))
                replacetext(newfile, 'MISALIGNNUM', '{}'.format(misalignnum))
                cmd = 'condor_submit ' + newfile
                os.system(cmd)