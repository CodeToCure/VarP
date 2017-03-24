#
# COPYRIGHT (C) 2017-2020 University of Maryland
#
"""
.. module:: plotGeneCoverage
    :platform: Unix, Windows, MacOSX
    :synopsis: Plots gene coverage  

.. moduleauthor:: Kunal Kundu (kkundu@umd.edu)

The module generates a gene coverage plot for a given gene

Module Dependency: 
 - Tabix
 - Numpy
 - Matplotlib
 - Varant

This module accepts only gVCF file.
It expects the gVCF should be Varant annotated.
"""


from gcn.lib.io.vcf import VCFParser
import os
import time
from tabix import Tabix
import matplotlib.pyplot as plt



def _load_capture(capture_file, qgene):
    '''Loads capture bed file'''
    cap_gene_coords = []
    s = open(capture_file)
    for line in s:
        line = line.strip()
        d = line.split('\t')
        if d[-1][:2] == 'rs':
            continue
        chrom, sp, ep, gene = d
        if gene != qgene:
            continue
        cap_gene_coords.append((chrom, int(sp), int(ep)))
    s.close()
    return cap_gene_coords


def _load_patient_capture(patient_capture_config, sample_ids):
    pc = {}
    s = open(patient_capture_config)
    for line in s:
        line = line.strip()
        if not line or line[0] in ['#', '-']:
            continue
        cver, pl = line.split('\t')
        pl = pl.split(',')
        for sample_id in sample_ids:
            if sample_id in pl:
                if cver not in pc:
                    pc[cver] = []
                pc[cver].append(sample_id)
    s.close()
    return pc


def _get_sample_ids(invcf):
    vcfo = VCFParser(invcf)
    return vcfo.samples
    vcfo.close()


def _get_coord_depth(invcf, capture_file_path, patient_capture_config, sample_ids, qgene):
    
    # Load the patient capture file
    pc = _load_patient_capture(patient_capture_config, sample_ids)
    
    print pc
    
    # Load the captured coordinates for each sample
    CAP = {'v2': '4-Hopkins_clinical_panel_capture_v2.bed',
           'v1b': '4-Hopkins_clinical_panel_capture_v1b.bed'}

    vcfo = Tabix(invcf)
    samples = _get_sample_ids(invcf)
    
    coord_depth = {}
    for cver, pl in pc.items():
        if cver == 'v2':
            continue
        capture_file = os.path.join(capture_file_path, CAP[cver])
        cap_exons = _load_capture(capture_file, qgene)
        ex_cnt = len(cap_exons)
        
        for idx, exon in enumerate(cap_exons):
            exno = idx + 1
            chrom, sp, ep = exon
            midpos = range(sp, ep)[(ep - sp)/2]
            for rec in vcfo.query(chrom, sp, ep):
                pos = int(rec[1])
                gt_info_raw = rec[9:]
                tsid = []
                for sid, gt_info in zip(samples, gt_info_raw):
                        if sid not in pl:
                            continue
                        tsid.append(sid)
                        if sid not in coord_depth:
                            coord_depth[sid] = [[], [], []]
                        if gt_info == './.':
                            dp = 0
                        else:
                            gt_d = gt_info.split(':')
                            if len(gt_d) <= 3:
                                dp = int(gt_d[-1])
                            else:
                                dp = int(gt_d[2])
                        coord_depth[sid][0].append(pos)
                        coord_depth[sid][1].append(dp)
                        if pos == midpos:
                            coord_depth[sid][2].append('CR-' + str(exno))
                        else:
                            coord_depth[sid][2].append('')
            
            if idx + 1 == ex_cnt:
                continue
            for sid in tsid:
                n = 105
                for fp in range(ep + 1, ep + 1 + n):
                    coord_depth[sid][0].append(fp)
                    coord_depth[sid][1].append(0)
                    coord_depth[sid][2].append('')
    return coord_depth

def gen_hex_colour_code():
    '''Generated hex color code'''
    
    import random
    colorhex = []
    for i in range(107):
        hx = ''.join([random.choice('0123456789ABCDEF') for x in range(6)])
        if hx not in colorhex:
            colorhex.append('#' + hx)
    return colorhex
      
def plot_gene_cov(coord_depth, qgene, colorhex):
    '''Plots the coverage'''
    
    temp = {}
    c = -1
    samples = coord_depth.keys()
    samples.sort()
    for sid in samples:
        v = coord_depth[sid]
        c += 1
        for a, b in zip(v[0], v[1]):
            if c == 0:
                temp[a] = []
            elif a not in temp:
                temp[a] = [0] * c
            temp[a].append(b)
            
        for k, e in temp.items():
            if len(e) != c + 1:
                temp[k].append(0)
    
    p = temp.keys()
    p.sort()
    
    yvals = []
    p = temp.keys()
    p.sort()
    for idx, k in enumerate(p):
        y = temp[k]
        if idx == 0:
            for e in range(len(y)):
                yvals.append([])
        for idx, yv in enumerate(y):
            yvals[idx].append(yv)
    
    
    xvals = range(1, len(p) + 1)
    
    for sid, yv, hx in zip(samples, yvals, colorhex):
        #print samples
        #print sid
        #hx = '#3F75A2' #Blue
        hx = '#936CA5'
        plt.plot(xvals, yv, label=sid, color=hx, lw=0.9, alpha=0.7)
        #plt.plot(xvals, yv, label=sid)
        #plt.legend()

    #plt.xlim(90694929, 90708737)
    #plt.ylim(0, 30)

    #plt.xticks(xvals, coord_depth[sid][2], rotation='vertical', fontsize=13)
    plt.ylabel('Read Depth', fontsize=18)
    plt.title('Coverage for %s gene' % qgene)
    plt.show()
        
        
def main(invcf, capture_file_path, genelist_config, patient_capture_config, sample_id, qgene, colorhex):
    coord_depth = _get_coord_depth(invcf, capture_file_path, patient_capture_config, sample_id, qgene)
    plot_gene_cov(coord_depth, qgene, colorhex)
    
    

if __name__ == '__main__':
    invcf = '/4-Hopkins_clinical_panel_SNV_dataset_v2.vcf.gz'
    capture_file_path = '/'
    genelist_config = '/genelist_config_v02.tsv'
    patient_capture_config = '/patient_capture_file_used.tsv'
    colorhex = gen_hex_colour_code()
    qgene = 'HYDIN'
    sample_ids = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P21', 'P22', 'P23', 'P24', 'P25', 'P26', 'P27', 'P28', 'P29', 'P30', 'P31', 'P32', 'P33', 'P34', 'P35', 'P36', 'P37', 'P38', 'P39', 'P40', 'P41', 'P42', 'P43', 'P44', 'P45', 'P46', 'P47', 'P48', 'P49', 'P50', 'P51', 'P52', 'P53', 'P54', 'P55', 'P56', 'P57', 'P58', 'P59', 'P60', 'P61', 'P62', 'P63', 'P64', 'P65', 'P66', 'P67', 'P68', 'P69', 'P70', 'P71', 'P72', 'P73', 'P74', 'P75', 'P76', 'P77', 'P78', 'P79', 'P80', 'P81', 'P82', 'P83', 'P84', 'P85', 'P86', 'P87', 'P88', 'P89', 'P90', 'P91', 'P92', 'P93', 'P94', 'P95', 'P96', 'P97', 'P98', 'P99', 'P100', 'P101', 'P102', 'P103', 'P104', 'P105', 'P106']
    t1 = time.time()
    main(invcf, capture_file_path, genelist_config, patient_capture_config, sample_ids, qgene, colorhex)
    t2 = time.time()
    print (t2 - t1) / 60,  'min'