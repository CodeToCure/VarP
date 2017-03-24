#
# COPYRIGHT (C) 2017-2020 University of Maryland
#
"""
.. module:: variantCallStat
    :platform: Unix, Windows, MacOSX
    :synopsis: Computes variant call related statistics 

.. moduleauthor:: Kunal Kundu (kkundu@umd.edu)

The module generates plots to illustrate - 
 * Ti/Tv ratio across samples
 * Homozygous/Heterozygous ratio across samples
 * Number of Common, Rare and Novel SNVs across samples
 * Number of Common, Rare and Novel Indels across samples
 * NUmber of no call sites and low quality sites across samples

Module Dependency: 
 - Tabix
 - Numpy
 - Matplotlib
 - Varant

This module accepts only gVCF file.
It expects the gVCF should be Varant annotated.
"""

from tabix import Tabix
from matplotlib import pyplot as plt
import numpy as np
import os
from gcn.lib.varann.vartype.varant.annotator import SNPAnnotation


def load_bed_file(capfile):
    '''Loads the given capture file'''
    capcoord = []
    s = open(capfile)
    for line in s:
        line = line.strip()
        if not line:
            continue
        d = line.split('\t')
        chrom, sp, ep, gene = d
        capcoord.append((chrom, int(sp), int(ep)))
    return capcoord


def _get_afrsamples(ethnicity_config):
    '''Returns a list of African samples'''
    afr_samples = []
    flag = False
    s = open(ethnicity_config)
    for line in s:
        line = line.strip()
        if not line:
            continue
        if flag == True:
            sid, e = line.split(':')
            if e.strip() == 'African':
                afr_samples.append(sid.strip())
        if line.startswith('Following'):
            flag = True
    return afr_samples


def _get_control_afrsamples(control_ethnicity_config):
    '''Returns a list of African samples'''
    control_afr_samples = []
    o = open(control_ethnicity_config)
    for line in o:
        line = line.strip()
        if not line:
            continue
        d = line.split('\t')
        if d[0] == 'sample':
            continue
        if d[2] == 'AFR':
            control_afr_samples.append(d[0])
    return control_afr_samples


def get_change_type(ref, alt):
    '''Given the Reference and Alternate allele,
    it return substitution type'''
    if len(ref) != len(alt):
        return 'NA'
    if ref in ['A', 'G'] and alt in ['A', 'G']:
        return 'ts'
    elif ref in ['C', 'T'] and alt in ['C', 'T']:
        return 'ts'
    else:
        return 'tv'


def get_vcfsamples(invcf):
    '''Given a VCF file it returns samples in the VCF'''
    for line in os.popen("zcat %s | head -2000 | grep '#CHROM'" % invcf):
        line = line.strip()
        h = line.split('\t')[:9]
        vcfsamples = line.split('\t')[9:]
        break
    return vcfsamples
    

def get_varcnt_gp(invcf, antvcf, capcoord, samples, varcnt, vao):
    '''Returns a dictionary where the keys are sampleid and values are a list
    of length 11 whose the elements are the count of - 
    * Nocall
    * VarLQ
    * GTLQ
    * REF allele
    * HetALT
    * HomALT
    * Ts
    * Tv
    * Common Variants (Based on 1000 Genomes)
    * Rare Variants (Based on 1000 Genomes)
    * Novel Variants (Based on 1000 Genomes)
    '''
    
    vcfsamples = get_vcfsamples(invcf)
    for sid in samples:
        varcnt[sid] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # NoCall, VarLQ, GTLQ, REF, HetALT, HomALT, Ts, Tv, Common, Rare, Novel
    vcfo = Tabix(invcf)
    antvcfo = Tabix(antvcf)
    for coord in capcoord:
        try:
            qdata = vcfo.query(*coord)
        except:
            qdata = None
        if qdata is None:
            continue
        for rec in qdata:
            filter = rec[6]
            gts = rec[9:]
            
            trans = vao.retrieve_genedef(rec[0], int(rec[1]), int(rec[1]))
            if not trans:  # indicating its an intergenic region
                continue
            
            for sid, gtinfo in zip(vcfsamples, gts):
                if sid in samples:
                    if gtinfo == '.' or './.' in gtinfo: # No Call
                        varcnt[sid][0] += 1
                    elif filter != 'PASS':
                        varcnt[sid][1] += 1 # VarLQ
                    elif '0/0' in gtinfo:
                        varcnt[sid][3] += 1 # REF
                    else:
                        gt, gq = gtinfo.split(':')[:2]
                        gq = int(gq)
                        a1, a2 = gt.split('/')
                        a1, a2 = int(a1), int(a2)
                        ref, alt = rec[3], rec[4].split(',')[a2-1]
                        ct = get_change_type(ref, alt)
                        arec = antvcfo.query(rec[0], int(rec[1])-1, int(rec[1]))
                        com, rare, novel = False, False, False
                        iflag = False
                        for e in arec:
                            if int(rec[1]) == int(e[1]):
                                kgaf = 0.0
                                if 'KGDB' in e[7]:
                                    kgaf = float(e[7].split('KGAF=')[1].split(';')[0].split(',')[a2 - 1])
                                exacaf = 0.0
                                if 'EXACDB' in e[7]:
                                    exacaf = float(e[7].split('EXACAF=')[1].split(';')[0].split(',')[a2 - 1])
                                    
                                if kgaf == 0.0 and exacaf == 0.0:
                                    novel = True
                                elif kgaf < 0.05 and exacaf < 0.05:
                                    rare = True
                                else:
                                    com = True
                                break
                        if gq < 30:   # GTLQ
                            varcnt[sid][2] += 1
                        elif a1 != a2:  # HetALT
                            varcnt[sid][4] += 1
                            if ct == 'ts':
                                varcnt[sid][6] += 1
                            elif ct == 'tv':
                                varcnt[sid][7] += 1
                            if com:
                                varcnt[sid][8] += 1
                            elif rare:
                                varcnt[sid][9] += 1
                            elif novel:
                                varcnt[sid][10] += 1
                                print sid, rec[0], rec[1], rec[3], rec[4]
                        elif a1 == a2:  # HomALT
                            varcnt[sid][5] += 1
                            if ct == 'ts':
                                varcnt[sid][6] += 1
                            elif ct == 'tv':
                                varcnt[sid][7] += 1
                            if com:
                                varcnt[sid][8] += 1
                            elif rare:
                                varcnt[sid][9] += 1
                            elif novel:
                                varcnt[sid][10] += 1
                                print sid, rec[0], rec[1], rec[3], rec[4]
    return varcnt


def plot_varcnt_by_type(varcnt, control_varcnt, samples_capv02, afr_samples, control_afr_samples, yl):
    '''Plot the variant call statistics'''
    
    control_commonlist = []
    control_rarelist_afr = []
    control_rarelist_nonafr = []
    
    for sid, val in control_varcnt.items():
        control_commonlist.append(val[4])
        if sid in control_afr_samples:
            control_rarelist_afr.append(val[5])
        else:
            control_rarelist_nonafr.append(val[5])
        
    legiden = {}
    for sid, val in varcnt.items():
        c, r, n = val[8], val[9], val[10]
        if sid in afr_samples:
            clr = '#1A5CEA'
            if sid in samples_capv02:
                m = '^'
                l = 'African;Capture v02'
            else:
                m = 'o'
                l = 'African;Capture v01'
        else:
            clr = '#EA831A'
            if sid in samples_capv02:
                m = '^'
                l = 'Non-African;Capture v02'
            else:
                m = 'o'
                l = 'Non-African;Capture v01'
        x1 = np.random.normal(0.5, 0.08, 1)
        p = plt.scatter(x1, c, s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
        if l not in legiden:
            legiden[l] = p
        x2 = np.random.normal(2, 0.08, 1)
        plt.scatter(x2, r, s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
        
        x3 = np.random.normal(3.8, 0.08, 1)
        plt.scatter(x3, n, s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
    
    plt.boxplot([control_commonlist, control_rarelist_afr, control_rarelist_nonafr], positions=[1, 2.5, 2.8])
    plt.xticks([0, 0.5, 1, 2, 2.5, 2.8, 3.8, 4.5], ['', 'HS', 'KGS','HS', 'KGS_AFR', 'KGS_NonAFR','HS', ''], rotation=50)
    plt.ylabel(yl)
    plt.ylim(0)
    t2 = []
    t1 = legiden.keys()
    t1.sort()
    for e in t1:
        t2.append(legiden[e])
    plt.grid(True)
    plt.show()


def plot_lq_nc_plot(varcnt, samples_capv02, afr_samples):
    '''Plot the No call sites and low quality variant count statistics'''
    
    legiden = {}
    for sid, val in varcnt.items():
        if sid in afr_samples:
            clr = '#1A5CEA'
            if sid in samples_capv02:
                m = '^'
                l = 'African;Capture v02'
            else:
                m = 'o'
                l = 'African;Capture v01'
        else:
            clr = '#EA831A'
            if sid in samples_capv02:
                m = '^'
                l = 'Non-African;Capture v02'
            else:
                m = 'o'
                l = 'Non-African;Capture v01'
        x1 = np.random.normal(val[0], 0.08, 1)
        p = plt.scatter(x1, val[1], s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
        if l not in legiden:
            legiden[l] = p
    
    plt.ylabel('# of No Call sites')
    plt.xlabel('# of low quality sites (not PASS in gVCF)')
    t2 = []
    t1 = legiden.keys()
    t1.sort()
    for e in t1:
        t2.append(legiden[e])
    plt.grid(True)
    plt.show()
    

def plot_ratio(varcnt, control_varcnt, samples_capv02, afr_samples, control_afr_samples):
    '''Plots the Ti/Tv ratio and HetALT/HomALT ratio across samples'''
    
    titv_ratio_list = []
    hethom_ratio_list = []
    sidlist = []
    control_titv_ratio_list = []
    control_hethom_ratio_list = []
    
    for s, val in varcnt.items():
        hethom_ratio = float(val[4]) / float(val[5])
        titv_ratio = float(val[6]) / float(val[7])
        hethom_ratio_list.append(hethom_ratio)
        titv_ratio_list.append(titv_ratio)
        sidlist.append(s)
    
    for s, val in control_varcnt.items():
        hethom_ratio = float(val[0]) / float(val[1])
        titv_ratio = float(val[2]) / float(val[3])
        control_hethom_ratio_list.append(hethom_ratio)
        control_titv_ratio_list.append(titv_ratio)
    
    legiden = {}
    for sid, t, h in zip(sidlist, titv_ratio_list, hethom_ratio_list):
        if sid in afr_samples:
            clr = '#1A5CEA'
            if sid in samples_capv02:
                m = '^'
                l = 'African;Capture v02'
            else:
                m = 'o'
                l = 'African;Capture v01'
        else:
            clr = '#EA831A'
            if sid in samples_capv02:
                m = '^'
                l = 'Non-African;Capture v02'
            else:
                m = 'o'
                l = 'Non-African;Capture v01'
        x1 = np.random.normal(0.5, 0.08, 1)
        p = plt.scatter(x1, t, s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
        if l not in legiden:
            legiden[l] = p
        x2 = np.random.normal(2, 0.08, 1)
        plt.scatter(x2, h, s=80, c=clr, alpha=0.8, marker=m, linewidths=0.5)
    
    plt.boxplot([control_titv_ratio_list, control_hethom_ratio_list], positions=[1, 2.5])
    plt.xticks([0, 0.5, 1, 2, 2.5, 3], ['', 'HS', 'KGS', 'HS', 'KGS', ''])
    plt.ylabel('Ratio')
    t2 = []
    t1 = legiden.keys()
    t1.sort()
    for e in t1:
        t2.append(legiden[e])
    plt.grid(True)
    plt.show()
    

def get_varcnt_control(controlvcf):
    '''Loads variants from the control vcf file'''
    
    control_varcnt = {}
    vcfsamples = get_vcfsamples(controlvcf)
    print '# of 1000 Genomes samples', len(vcfsamples)
    for sid in vcfsamples:
        control_varcnt[sid] = [0, 0, 0, 0, 0, 0, 0]  # HET, HOM, Ts, Tv, Common, Rare, Novel
    s = os.popen("zcat %s" % controlvcf)
    for rec in s:
        if rec[0] == '#' :
            continue
        rec = rec.strip()
        rec = rec.split('\t')
        for sid, gtinfo in zip(vcfsamples, rec[9:]):
            if '.' in gtinfo or gtinfo in ['0/0', '0|0', '0']:
                continue
            gt = gtinfo
            a1, a2 = gt.split('|')
            a1, a2 = int(a1), int(a2)
            ref, alt = rec[3], rec[4].split(',')[a2-1]
            ct = get_change_type(ref, alt)
            com, rare = False, False
            af = float(rec[7].split(';AF=')[1].split(';')[0].split(',')[a2 - 1])
            #print af
            if af >= 0.05:
                com = True
            #elif af 
            else:
                rare = True  
            if a1 != a2:  # HetALT
                control_varcnt[sid][0] += 1
                if ct == 'ts':
                    control_varcnt[sid][2] += 1
                elif ct == 'tv':
                    control_varcnt[sid][3] += 1
                if com:
                    control_varcnt[sid][4] += 1
                elif rare:
                    control_varcnt[sid][5] += 1
            elif a1 == a2:  # HomALT
                control_varcnt[sid][1] += 1
                if ct == 'ts':
                    control_varcnt[sid][2] += 1
                elif ct == 'tv':
                    control_varcnt[sid][3] += 1
                if com:
                    control_varcnt[sid][4] += 1
                elif rare:
                    control_varcnt[sid][5] += 1
    return control_varcnt
                
        
    

def main(snvvcf, indelvcf, snvantvcf, snvcontrolvcf, indelcontrolvcf, capv01, capv02, samples_capv01, samples_capv02, ethnicity_config, control_ethnicity_config):
    
    vao = SNPAnnotation('REFGENE')
 
    cap1coord = load_bed_file(capv01)
    cap2coord = load_bed_file(capv02)
    afr_samples = _get_afrsamples(ethnicity_config)
    control_afr_samples = _get_control_afrsamples(control_ethnicity_config)
    
    ## FOR SNV ##
    # Load the SNV Count data for all captures
    snvvarcnt = {}
    snvvarcnt = get_varcnt_gp(snvvcf, snvantvcf, cap1coord, samples_capv01, snvvarcnt, vao)
    snvvarcnt = get_varcnt_gp(snvvcf, snvantvcf, cap2coord, samples_capv02, snvvarcnt, vao)
    
    # Load the SNV variant Count data for Control
    snv_control_varcnt = get_varcnt_control(snvcontrolvcf)
    
    
    ## FOR INDEL ##
    # Load the INDEL Count data for all captures
    antvcf = indelvcf
    indelvarcnt = {}
    indelvarcnt = get_varcnt_gp(indelvcf, antvcf, cap1coord, samples_capv01, indelvarcnt, vao)
    indelvarcnt = get_varcnt_gp(indelvcf, antvcf, cap2coord, samples_capv02, indelvarcnt, vao)
    
    # Load the INDEL variant Count data for Control
    indel_control_varcnt = get_varcnt_control(indelcontrolvcf)
    
    
    
    # Display the Loaded Variant Count data in csv format
    print '******'
    for s, val in indelvarcnt.items():
        sum(val[:6])
        if val[6] == 0 or val[7] == 0:
            titv_ratio = '0.0'
        else:
            titv_ratio = str(float(val[6]) / float(val[7]))
        print s, '-', ','.join([str(e) for e in val] + [titv_ratio])
    print '******'
    
    # Display the Loaded Variant Count data in csv format
    print '########'
    for s, val in snvvarcnt.items():
        sum(val[:6])
        if val[6] == 0 or val[7] == 0:
            titv_ratio = '0.0'
        else:
            titv_ratio = str(float(val[6]) / float(val[7]))
        print s, '-', ','.join([str(e) for e in val] + [titv_ratio])
    print '########'
    
    # Plot No call versus Low Call sites
    plot_lq_nc_plot(snvvarcnt, samples_capv02, afr_samples)
    
    #Plot the INDEL - Common, Rare, Novel variants
    yl = '# of Indels'
    plot_varcnt_by_type(indelvarcnt, indel_control_varcnt, samples_capv02, afr_samples, control_afr_samples, yl)
    
    #Plot the SNV - Common, Rare, Novel variants
    yl = '# of SNVs'
    plot_varcnt_by_type(snvvarcnt, snv_control_varcnt, samples_capv02, afr_samples, control_afr_samples, yl)
    
    # Plot Ts/Tv ratio and HET/HOM ratio
    plot_ratio(snvvarcnt, snv_control_varcnt, samples_capv02, afr_samples, control_afr_samples)
    
                            
if __name__ == '__main__':
    snvvcf = '/4-Hopkins_clinical_panel_SNV_dataset_v2.vcf.gz'
    indelvcf = '/4-Hopkins_clinical_panel_InDel_varant_refgene.vcf.gz'
    snvcontrolvcf = '/control_variants_all.vcf.gz'
    indelcontrolvcf = '/indel_control_variants_all.vcf.gz'
    snvantvcf = '/4-Hopkins_clinical_panel_SNV_varant_refgene.vcf.gz'
    indelantvcf = '/4-Hopkins_clinical_panel_InDel_varant_refgene.vcf.gz'
    capv01 = '/4-Hopkins_clinical_panel_capture_v1b.bed'
    capv02 = '/4-Hopkins_clinical_panel_capture_v2paper.bed'
    samples_capv01 = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P13',
                      'P14', 'P15', 'P16', 'P17', 'P18', 'P19', 'P20', 'P23', 'P24', 'P25', 'P26',
                      'P28', 'P29', 'P31', 'P32', 'P33', 'P35', 'P36', 'P37', 'P38', 'P40',
                      'P41', 'P42', 'P43', 'P44', 'P45', 'P46', 'P48', 'P49', 'P51',
                      'P52', 'P53', 'P54', 'P56', 'P57', 'P58', 'P59', 'P60', 'P61', 'P62',
                      'P63', 'P64', 'P65', 'P66', 'P67', 'P68', 'P69', 'P70', 'P71', 'P72', 'P73', 
                      'P74', 'P75', 'P76', 'P77', 'P78', 'P79', 'P80', 'P81', 'P82', 'P83', 'P84',
                      'P85', 'P86', 'P87', 'P88', 'P89', 'P90', 'P91', 'P92', 'P93', 'P94', 'P95', 
                      'P96', 'P97', 'P98', 'P99', 'P100', 'P101', 'P102', 'P103', 'P104', 'P105', 'P106']
    samples_capv02 = ['P12', 'P21', 'P22', 'P27', 'P30', 'P34', 'P39', 'P47', 'P50', 'P55']
    ethnicity_config = '/4-Hopkins_clinical_panel_SNV_predicted_ethnicity_KK.tsv'
    control_ethnicity_config = '/integrated_call_samples_v3.20130502.ALL.panel'
    control_ped_file = '/20130606_g1k.ped'
    main(snvvcf, indelvcf, snvantvcf, snvcontrolvcf, indelcontrolvcf, capv01, capv02, samples_capv01, samples_capv02, ethnicity_config, control_ethnicity_config)


