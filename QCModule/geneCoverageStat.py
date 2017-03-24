#
# COPYRIGHT (C) 2017-2020 University of Maryland
#
"""
.. module:: geneCoverageStat
    :platform: Unix, Windows, MacOSX
    :synopsis: Computes gene coverage statistics 

.. moduleauthor:: Kunal Kundu (kkundu@umd.edu)

The module computed average read depth per gene across samples.
It also identifies region with anomalous coverage.

Module Dependency: 
 - Tabix
 - Numpy
 - Varant

This module accepts only gVCF file.
It expects the gVCF should be Varant annotated.
"""

from gcn.lib.io.vcf import VCFParser
import numpy as np
import os
import time
from tabix import Tabix


def _load_capture(capture_file, gl):
    '''Loaded the Capture bed file''' 
    
    cap_gene_coords = {}
    s = open(capture_file)
    for line in s:
        line = line.strip()
        d = line.split('\t')
        if d[-1][:2] == 'rs':
            continue
        chrom, sp, ep, gene = d
        if gene not in cap_gene_coords:
            cap_gene_coords[gene] = []
        cap_gene_coords[gene].append((chrom, int(sp), int(ep)))
    s.close()
    return cap_gene_coords


def _load_genelist(genelist_config):
    gl = []
    s = open(genelist_config)
    for line in s:
        line = line.strip()
        if not line or line[0] in ['#', '-']:
            continue
        gene = line.split('\t')[0]
        gl.append(gene)
    s.close()
    gl = list(set(gl))
    gl.sort()
    return gl


def _load_patient_capture(patient_capture_config):
    pc = {}
    s = open(patient_capture_config)
    for line in s:
        line = line.strip()
        if not line or line[0] in ['#', '-']:
            continue
        cver, pl = line.split('\t')
        pl = pl.split(',')
        pc[cver] = pl
    s.close()
    return pc


def _get_sample_ids(invcf):
    '''Return sample ids present in the vcf file'''
    vcfo = VCFParser(invcf)
    return vcfo.samples
    vcfo.close()


def _get_coord_depth(invcf, capture_file_path, genelist_config, patient_capture_config):
    #global capcoord_dp
    
    # Load the genelist
    gl = _load_genelist(genelist_config)
    #gl = ['CCDC39', 'CSF2RA', 'CTC1', 'DNM1L', 'FBN1', 'HSD17B4', 'HYDIN', 'TERT']
    
    # Load the patient capture file
    pc = _load_patient_capture(patient_capture_config)
    
    # Load the captured coordinates for each sample
    CAP = {'v2': '4-Hopkins_clinical_panel_capture_v2paper.bed',
           'v1b': '4-Hopkins_clinical_panel_capture_v1b.bed'}

    cap_exon_mean = {}
    vcfo = Tabix(invcf)
    samples = _get_sample_ids(invcf)
    for cver, pl in pc.items():
        capture_file = os.path.join(capture_file_path, CAP[cver])
        cap_gene_coords = _load_capture(capture_file, gl)
        for gene, cap_exons in cap_gene_coords.items():
            for idx, exon in enumerate(cap_exons):
                exno = idx + 1
                chrom, sp, ep = exon
                key = (chrom, sp, ep, exno)
                temp = {}
                for rec in vcfo.query(chrom, sp, ep):
                    gt_info_raw = rec[9:]
                    for sid, gt_info in zip(samples, gt_info_raw):
                        if sid not in pl:
                            continue
                        if sid not in temp:
                            temp[sid] = []
                        if gt_info == './.':
                            temp[sid].append(0)
                        else:
                            gt_d = gt_info.split(':')
                            if len(gt_d) <= 3:
                                temp[sid].append(int(gt_d[-1]))
                            else:
                                temp[sid].append(int(gt_d[2]))
                # Compute the mean
                for sid, val in temp.items():
                    m = np.average(val)
                    if sid not in cap_exon_mean:
                        cap_exon_mean[sid] = {}
                    if gene not in cap_exon_mean[sid]:
                        cap_exon_mean[sid][gene] = {}
                    cap_exon_mean[sid][gene][key] = m
    return cap_exon_mean
        
      
        
        
        
def main(invcf, capture_file_path, genelist_config, patient_capture_config):
    n = 1
    om = open(os.path.join(capture_file_path, 'output_paper.mean.tsv'), 'w')
    od = open(os.path.join(capture_file_path, 'output_paper.sd.tsv'), 'w')
    oexdup = open(os.path.join(capture_file_path, 'paper_ex_dup_' + str(n) + '.tsv'), 'w')
    oexdel = open(os.path.join(capture_file_path, 'paper_ex_del_' + str(n) + '.tsv'), 'w')
    
    np.nan_to_num(0)
    cap_exon_mean = _get_coord_depth(invcf, capture_file_path, genelist_config, patient_capture_config)

    gl = _load_genelist(genelist_config)
    #gl += ['DHODH', 'TRIM37', 'EFTUD2', 'AMACR', 'CAT', 'AGXT']
    gl.sort()
    sample_ids = cap_exon_mean.keys()
    sample_ids.sort()
    d = {}
    m = {}
    ex_del = {}
    ex_dup = {}
    for idx, sid in enumerate(sample_ids):
        for g in gl:
            if g not in d:
                d[g] = []
                m[g] = []
            
            if g not in cap_exon_mean[sid]:
                d[g].append('NA')
                m[g].append('NA')
                continue
            
            y = cap_exon_mean[sid][g].values()
            sd = np.std(y)
            meancov_gene = np.average(y)
            meancov_gene = np.nan_to_num(meancov_gene)
            #print sid, g, meancov_gene, sd
            nstd = n * sd
            d[g].append(sd)
            m[g].append(meancov_gene)
            for key, ad in cap_exon_mean[sid][g].items():
                if ad > meancov_gene + nstd:
                    chrom, sp, ep, exno = key
                    k = '_'.join([g, chrom, str(sp), str(ep), 'Exon-' + str(exno)])
                    if idx == 0:
                        ex_dup[k] = []
                    elif k not in ex_dup:
                        ex_dup[k] = [0] * idx
                    ex_dup[k].append(1)
                elif ad < meancov_gene - nstd:
                    chrom, sp, ep, exno = key
                    k = '_'.join([g, chrom, str(sp), str(ep), 'Exon-' + str(exno)])
                    if idx == 0:
                        ex_del[k] = []
                    elif k not in ex_del:
                        ex_del[k] = [0] * idx
                    ex_del[k].append(1)
        for a, b in ex_dup.items():
            if len(b) != idx + 1:
                ex_dup[a].append(0)
        for a, b in ex_del.items():
            if len(b) != idx + 1:
                ex_del[a].append(0)
                
    # Print the matrix to file
    p = '\t'.join([''] + sample_ids) + '\n'
    om.write(p)
    od.write(p)
    oexdup.write(p)
    oexdel.write(p)
    for x, y in m.items():
        p = '\t'.join([str(e) for e in [x] + y]) + '\n'
        om.write(p)
    
    for x, y in d.items():
        p = '\t'.join([str(e) for e in [x] + y]) + '\n'
        od.write(p)
    
    kl = ex_dup.keys()
    kl.sort()
    for k in kl:
        p = '\t'.join([str(e) for e in [k] + ex_dup[k]]) + '\n'
        oexdup.write(p)
        
    kl = ex_del.keys()
    kl.sort()
    for k in kl:
        p = '\t'.join([str(e) for e in [k] + ex_del[k]]) + '\n'
        oexdel.write(p)
        
    om.close()
    od.close()
    oexdup.close()
    oexdel.close()

if __name__ == '__main__':
    invcf = '/4-Hopkins_clinical_panel_SNV_dataset_v2.vcf.gz'
    capture_file_path = '/'
    genelist_config = '/genelist_config_vpaper.tsv'
    patient_capture_config = '/patient_capture_file_used.tsv'
    t1 = time.time()
    main(invcf, capture_file_path, genelist_config, patient_capture_config)
    t2 = time.time()
    print (t2 - t1) / 60,  'min'