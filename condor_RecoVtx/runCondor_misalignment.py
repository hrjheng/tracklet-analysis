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

dphicut_bin = 3
dzcut_bin = 10
gap_north = 3.5
gap_atClips = 1.5
gap_upper_min = gap_atClips / 2.
gap_upper_max = gap_north - (gap_atClips / 2.)
Nbinscan_gap = 100
centshift_min = -0.5
centshift_max = 0.5
Nbinscan_centshift = 50
gap_upper = [gap_upper_min + (gap_upper_max - gap_upper_min) / Nbinscan_gap * i for i in range(0, Nbinscan_gap+1)]
cent_shift = [centshift_min + (centshift_max - centshift_min) / Nbinscan_centshift * i for i in range(0, Nbinscan_centshift+1)]

for i, gapupper in enumerate(gap_upper):
    for j, centshift in enumerate(cent_shift):
        print('(dPhiCut,dZCut,gap_north, gap_upper, cent_shift)=({:.2f},{:.2f},{:.2f},{:.2f},{:.2f})'.format(dphicut_bin*0.01, dzcut_bin*0.01, gap_north, gapupper, centshift))
        gapupper_str = '{:.2f}'.format(gapupper).replace('.', 'p')
        centshift_str = '{:.2f}'.format(centshift).replace('.', 'p')
        newfile = 'condor_dphicut{}_dzcut{}_gapnorth3p5_gapupperbin{}_centshiftbin{}.job'.format(dphicut_bin, dzcut_bin, i, j)
        os.system('cp condor.job {}'.format(newfile))
        replacetext(newfile, 'DPHICUT', '{}'.format(dphicut_bin))
        replacetext(newfile, 'DZCUT', '{}'.format(dzcut_bin))
        replacetext(newfile, 'GAPNORTH', '{:.1f}'.format(gap_north))
        replacetext(newfile, 'GAPUPPER', '{:.2f}'.format(gapupper))
        replacetext(newfile, 'CENTSHIFT', '{:.2f}'.format(centshift))

        cmd = 'condor_submit ' + newfile
        os.system(cmd)


list_gap_north = [1.5 + i * 0.1 for i in range(21)]
for i, gapnorth in enumerate(list_gap_north):
    gapupper = gapnorth / 2.
    print('(dPhiCut, dZCut, gap_north, gap_upper, cent_shift)=({:.2f},{:.2f},{:.2f},{:.2f},{:.2f})'.format(dphicut_bin*0.01, dzcut_bin*0.01, gapnorth, gapnorth / 2., 0))
    gapnorth_str = '{:.1f}'.format(gapnorth).replace('.', 'p')
    gapupper_str = '{:.2f}'.format(gapupper).replace('.', 'p')
    newfile = 'condor_dphicut{}_dzcut{}_gapnorth{}_gapupperbin{}_centshift0.job'.format(dphicut_bin, dzcut_bin, gapnorth_str, gapupper_str)
    os.system('cp condor.job {}'.format(newfile))
    replacetext(newfile, 'DPHICUT', '{}'.format(dphicut_bin))
    replacetext(newfile, 'DZCUT', '{}'.format(dzcut_bin))
    replacetext(newfile, 'GAPNORTH', '{:.1f}'.format(gapnorth))
    replacetext(newfile, 'GAPUPPER', '{:.2f}'.format(gapupper))
    replacetext(newfile, 'CENTSHIFT', '{:.1f}'.format(0))

    cmd = 'condor_submit ' + newfile
    os.system(cmd)