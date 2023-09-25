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

for i in range(0, 30):
    for j in range(0, 30):
        dphicut = 0.01 * (i+1)
        dzcut = 0.01 * (j+1)
        print('(dPhiCut,dZCut)=({:.2f},{:.2f})'.format(dphicut, dzcut))
        newfile = 'condor_dphicut{}_dzcut{}.job'.format(i+1,j+1)
        os.system('cp condor.job {}'.format(newfile))
        replacetext(newfile, 'DPHICUT', '{}'.format(i+1))
        replacetext(newfile, 'DZCUT', '{}'.format(j+1))
        cmd = 'condor_submit ' + newfile
        os.system(cmd)