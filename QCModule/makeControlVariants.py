#
# COPYRIGHT (C) 2017-2020 University of Maryland
#
"""
.. module:: makeControlVariants
    :platform: Unix, Windows, MacOSX
    :synopsis: Extraction of control variants from 1000 Genomes Phase-3 data

.. moduleauthor:: Kunal Kundu (kkundu@umd.edu)

The module generates a TSV file comprised of variants reported in 1000 Genomes
Project Phase 3 data set in given genomic regions. The genomic regions are provided in 
the form of a bed file.

Module Dependency: 
 - Tabix
 
"""
import os
import random
from tabix import Tabix


def _get_random_samples_info(outdir, num_afr_smpls, num_nonafr_smpls):
    """Internal Method. Gets a list of randomized African and 
    randomized non-African samples present in 1000 Genomes Project 
    Phase 3 dataset. 

    Args:
        outdir (str): Path where 1000 Genomes file will be downloaded
        num_afr_smpls (str): Number of random African samples
        num_nonafr_smpls (str): Number of random nonAfrican samples

    Returns:
        sample_info (dictionary): Keys are 'AFR' whose value is a list 
                                  of African samples and 'nAFR' whose 
                                  value is a list of non-African samples
    """
    samples_info = {'AFR': [], 'nAFR': []}
    ofp = os.path.join(outdir, 'integrated_call_samples_v3.20130502.ALL.panel')
    if not os.path.exists(ofp):
        fp = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
        cmd = "wget %s  -O %s" % (fp, ofp)
        os.system(cmd)
    o = open(ofp)
    for line in o:
        line = line.strip()
        if not line:
            continue
        d = line.split('\t')
        if d[0] == 'sample':
            continue
        if d[2] == 'AFR':
            samples_info['AFR'].append(d[0])
        else:
            samples_info['nAFR'].append(d[0])
    
    samples_info['AFR'] = random.sample(samples_info['AFR'], num_afr_smpls)
    samples_info['nAFR'] = random.sample(samples_info['nAFR'], num_nonafr_smpls)
    return  samples_info


def _get_samples_info(outdir):
    """Internal Method. Gets a list of African and non-African samples 
    present in 1000 Genomes Project Phase 3 dataset. 

    Args:
        outdir (str): Path where 1000 Genomes file will be downloaded

    Returns:
        sample_info (dictionary): Keys are 'AFR' whose value is a list 
                                  of African samples and 'nAFR' whose 
                                  value is a list of non-African samples
    """
    samples_info = {'AFR': [], 'nAFR': []}
    ofp = os.path.join(outdir, 'integrated_call_samples_v3.20130502.ALL.panel')
    if not os.path.exists(ofp):
        fp = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
        cmd = "wget %s  -O %s" % (fp, ofp)
        os.system(cmd)
    o = open(ofp)
    for line in o:
        line = line.strip()
        if not line:
            continue
        d = line.split('\t')
        if d[0] == 'sample':
            continue
        if d[2] == 'AFR':
            samples_info['AFR'].append(d[0])
        else:
            samples_info['nAFR'].append(d[0])
    #os.system('rm -rf %s' % ofp)
    return samples_info


def _generate_TSV(outdir, samples, capcoord):
    """Internal Method. Generates two tsv files (one for snv and 
    other for indel) that contains the variants reported in 1000
    Genomes Project Phase 3 dataset for given genomic regions. 

    Args:
        outdir (str): Path where TSV files will be created
        samples (list): List of samples for which the variants
                        will be extracted from 1000 Genomes dataset.
        capcoord (list): List of bed format coordinates wrapped as 
                         tuples. For eg. 
                         [(chrom, start_pos, end_pos), ...]

    Returns:
        None
    """
    out_snpfile = os.path.join(outdir, 'snp_control_variants_all.tsv')
    outsnp = open(out_snpfile, 'w')

    out_indelfile = os.path.join(outdir, 'indel_control_variants_all.tsv')
    outindel = open(out_indelfile, 'w')
    
    chroms = [str(e) for e in range(1, 23)] + ['X', 'Y']
    hflag = False
    for chrom in chroms:
        if chrom == 'X':
            fp1 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
            ofp1 = os.path.join(outdir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')
            fp2 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi'
            ofp2 = os.path.join(outdir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi')
        else:
            fp1 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % chrom
            ofp1 = os.path.join(outdir, 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' % chrom)
            fp2 = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi' % chrom
            ofp2 = os.path.join(outdir, 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi' % chrom)
        
        if not os.path.exists(ofp1):
            cmd = "wget %s -O %s" % (fp1, ofp1)
            print cmd
            os.system(cmd)
        if not os.path.exists(ofp2):
            cmd = "wget %s -O %s" % (fp2, ofp2)
            print cmd
            os.system(cmd)
        
        for line in os.popen("zcat %s | head -260 | grep '#CHROM'" % ofp1):
            line = line.strip()
            h = line.split('\t')[:9]
            kgsamples = line.split('\t')[9:]
        
        if hflag == False:
            h += samples
            outsnp.write('\t'.join(h) + '\n')
            outindel.write('\t'.join(h) + '\n')
            hflag = True
        
        vcfo = Tabix(ofp1)
        for coord in capcoord:
            if coord[0] != chrom:
                continue
            for rec in vcfo.query(*coord):
                if '<' in rec[4] or ',' in rec[4]:
                    continue
                
                if len(rec[3]) == len(rec[4]):  # to make Indel control set
                    out = outsnp
                
                if len(rec[3]) != len(rec[4]):  # to make SNV control set
                    out = outindel
            
                sv = []
                d = dict(zip(kgsamples, rec[9:]))
                temp = set([])
                for sid in samples:
                    sv.append(d[sid])
                    temp.add(d[sid])
                if len(temp) == 1 and list(temp)[0] in ['0|0', '0/0']:
                    continue
                out.write('\t'.join(rec[:9] + sv) + '\n')
        #os.system('rm -rf %s' % ofp1)
        #os.system('rm -rf %s' % ofp2)
                
        
def _load_capture_file(capture_file):
    """Internal Method. Loads the capture bed file. 

    Args:
        capture_file (str): Path to capture bed file

    Returns:
        capcoord (list): List of bed format coordinates wrapped as 
                         tuples. For eg. 
                         [(chrom, start_pos, end_pos), ...]
    """
    capcoord = []
    s = open(capture_file)
    for line in s:
        line = line.strip()
        if not line:
            continue
        d = line.split('\t')
        chrom, sp, ep, gene = d
        if gene[2:].isdigit():
            continue
        capcoord.append((chrom, int(sp), int(ep)))
    return capcoord


def main(outdir, capture_file, random_flag, num_afr, num_nonafr):
    """Main method to generate the TSV files  

    Args:
        outdir (str): Path where TSV files will be created
        capture_file (str): Path to capture bed file
        random_flag (boolean): Flag to denote of the 1000 Genomes 
                               samples will be randomized for extracting 
                               variants
        num_afr (integer): Number of African samples randomly selected
        num_nonafr(integer): Number of Non-Afican samples randomly selected

    Returns:
        None
    """
    if random_flag == True:
        samples_info = _get_random_samples_info(outdir, num_afr, num_nonafr)
    else:
        samples_info = _get_samples_info(outdir)
    samples = samples_info['AFR'] + samples_info['nAFR']
    samples.sort()
    capcoord = _load_capture_file(capture_file)
    _generate_TSV(outdir, samples, capcoord)


if __name__ == '__main__':
    import argparse
    desc = 'Script to generate TSV files (snp and indel) containg \
    1000 Genomes variants present in the coordinates of the given bed file'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-o', '--outputDir', dest='outdir', type=str,
                        help='Folder path where TSV files will be created')
    parser.add_argument('-c', '--capturebed', dest='capture_file', type=str,
                        help='Path to Capture bed file')
    parser.add_argument('-r', '--random', dest='random_flag', type=str,
                        default=False, help='Flag to denote of the 1000 Genomes\
                        samples will be randomized for extracting variants. Enter: \
                        True or False. Default: False')
    parser.add_argument('-a', '--numafr', dest='num_afr', type=str,
                        help='If random is True, then specify \
                        number of african samples you want to consider')
    parser.add_argument('-n', '--numnonafr', dest='num_nonafr', type=str,
                        help='If random is True, then specify \
                        number of non-african samples you want to consider')
    args = parser.parse_args()
    if not args.outdir:
        print "Pls enter the 'outdir' value"
        exit(1)
    if not args.capture_file:
        print "Pls enter the 'capturebed' value"
        exit(1)
    if args.random_flag == 'True':
        rflag = True
        if not args.num_afr:
            print "Pls enter the 'num_afr' value"
            exit(1)
        else:
            num_afr = args.num_afr
        if not args.num_nonafr:
            print "Pls enter the 'num_nonafr' value"
            exit(1)
        else:
            num_nonafr = args.num_nonafr
    elif args.random_flag == 'False':
        rflag = False
        num_afr, num_nonafr = 0, 0
    else:
        rflag = False
        num_afr, num_nonafr = 0, 0
    main(args.outdir, args.capture_file, rflag, num_afr, num_nonafr)