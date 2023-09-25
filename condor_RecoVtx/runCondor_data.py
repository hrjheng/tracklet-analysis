import time
import os
import sys
import re
from argparse import ArgumentParser
from datetime import datetime


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

# Setup
dphicut_bin = 3
dzcut_bin = 10
gapnorth = gap_atClips = 1.5
gapupper = gap_atClips / 2.
centshift = 0

print('Run on (mock)data: (dPhiCut, dZCut, gapnorth, gap_upper, cent_shift)=({:.2f},{:.2f},{:.2f},{:.2f},{:.2f})'.format(dphicut_bin*0.01, dzcut_bin*0.01, gapnorth, gapupper, centshift))
gapnorth_str = '{:.2f}'.format(gapnorth).replace('.', 'p')
gapupper_str = '{:.2f}'.format(gapupper).replace('.', 'p')
centshift_str = '{:.2f}'.format(centshift).replace('.', 'p')
newfile = 'condor_rundata_dphicut{}_dzcut{}_gapnorth{}_gapupper{}_centshift{}.job'.format(dphicut_bin, dzcut_bin, gapnorth_str, gapupper_str, centshift_str)
os.system('cp condor.job {}'.format(newfile))
replacetext(newfile, 'ISDATA', '{}'.format(1))
replacetext(newfile, 'DPHICUT', '{}'.format(dphicut_bin))
replacetext(newfile, 'DZCUT', '{}'.format(dzcut_bin))
replacetext(newfile, 'GAPNORTH', '{:.1f}'.format(gapnorth))
replacetext(newfile, 'GAPUPPER', '{:.2f}'.format(gapupper))
replacetext(newfile, 'CENTSHIFT', '{:.2f}'.format(centshift))

cmd = 'condor_submit ' + newfile
os.system(cmd)