#
# COPYRIGHT (C) 2017-2020 University of Maryland
#
"""
.. module:: predictDisease
    :platform: Unix, Windows, MacOSX
    :synopsis: Computes variant call related statistics 

.. moduleauthor:: Kunal Kundu (kkundu@umd.edu)

The module priotizes variants per sample and further predicts the disease class

Module Dependency: 
 - Tabix
 - Varant

This module accepts only VCF file.
It expects the gVCF should be Varant annotated.
"""


from gcn.lib.io.vcf import VCFParser
from gcn.lib.varann.vartype.varant import varant_parser as vp
from gcn.lib.databases.genedef import GeneDef
from gcn.data.complement import COMPLEMENT
from tabix import Tabix


CLNSIG_MAP = {0: 'Uncertain significance',
              1: 'not provided',
              2: 'Benign', 
              3: 'Likely benign',
              4: 'Likely pathogenic',
              5: 'Pathogenic',
              6: 'drug response',
              7: 'histocompatibility', 
              255: 'other'}



def isRare(altid, info, thres_af):
    '''Checks if the variant is rare'''
    
    thres = thres_af
    af = 0.00
    flag = False
    '''if 'KGAF' in info:
        af = float(info['KGAF'][altid - 1])
    else:
        af = 0.0'''
    
    if af == 0.0:
        if 'EXACAF' in info:
            af = float(info['EXACAF'][altid - 1])
    '''if af == 0.0:
        if 'ESPAF' in info:
            af = float(info['ESPAF'][altid - 1])'''
    
    if af <= thres:
        flag = True
    else:
        flag = False
    return (af, flag)
    
    
def _genotype_check(gt, gq):
    '''Return False if the GQ < 30'''
    
    if '.' in gt or gt == '0/0':
        return False
    elif float(gq) < 30:
        return False
        #return True
    else:
        return True


def _is_HQVar(filter_ant):
    '''Checks if variant is PASS'''
    
    if 'PASS' in filter_ant:
        return True
    else:
        return False
    #return True


def _is_eQTLs(altid, info):
    '''Checks if variant is an eQTL'''
    
    if 'RegulomeScore' in info:
        if info['RegulomeScore'][altid - 1][0] == '1':
            return True
        else:
            return False
    else:
        return False

def _is_PASnv(ta):
    if 'NonSyn' in ta.mutation or 'StartLoss' in ta.mutation:
        return True
    else:
        return False

def _is_NonSense(ta):
    if 'Stop' in ta.mutation:
        return True
    else:
        return False

def _is_Splicing(ta):
    if ta.splice in ['SpliceDonor', 'SpliceAcceptor']:
    #if 'Splice' in ta.splice:
        return True
    else:
        return False


def _is_PAIndel(ta):
    if 'Frame' in ta.mutation:
        return True
    else:
        return False

def _is_Intronic(ta):
    if 'CodingIntronic' in ta.region:
        return True
    else:
        return False
    
def _is_UTR(ta):
    if 'UTR' in ta.region:
        return True
    else:
        return False

def _is_Damaging(altid, info, ta, snps3d_pred, nmethod):
    profile_pred, stability_pred = snps3d_pred[1], snps3d_pred[3]
    
    if '_' in ta.sift:
        s = ta.sift.split('_')[0]
    else:
        s = ''
    if '_' in ta.pp2:
        pp2 = ta.pp2.split('_')[0]
    else:
        pp2 = ''
    cadd = float(info['CADD_phred'][altid - 1])
    
    if nmethod == 1:
        if s == 'D' or pp2 in ['PP2D', 'PP2PD'] or cadd >= 15.0 or \
        'High' in profile_pred or 'High' in stability_pred:
            return True
        else:
            return False
    elif nmethod == 2:
        if (s == 'D' and pp2 in ['PP2D', 'PP2PD']) or \
        (s == 'D' and cadd >= 15.0) or (pp2 in ['PP2D', 'PP2PD'] \
                                        and cadd >= 15.0):
            return True
        else:
            return False
    elif nmethod == 3:
        if s == 'D' and pp2 in ['PP2D', 'PP2PD'] and cadd >= 15.0:
            return True
        else:
            return False
    else:
        return False

    
def get_dbscSNV_ant(chrom, pos, ref, alt):
    '''Checks if the variant affects splicing'''
    f = '/dbscSNV/dbscSNV1.1.chr' + chrom + '.gz'
    fo = Tabix(f)
    for rec in fo.query(chrom, pos - 1, pos + 1):
        if int(rec[1]) == pos and rec[2] == ref and rec[3] == alt:
            ada_score, rf_score = rec[-2], rec[-1].strip()
            if rf_score and rf_score != '.':
                rf_score = ''
            return float(ada_score), rf_score
    return '', ''
        
        

def _load(invcf, thres_af, nmethod, data = None, sc = 1):
    if not data:
        data = {}
    vcfs = VCFParser(invcf)
    samples = vcfs.samples
    for rec in vcfs:
        vcfs.parseinfo(rec)
        vcfs.parsegenotypes(rec)
        
        if not _is_HQVar(rec.filter):  # Checks variant is PASS
            continue
        
        for sid in samples:
            if sid not in data:
                data[sid] = {}
            
            gi = rec[sid]
            gt, gq = gi.GT, gi.GQ
            
            # Checks if genotype is not reference or GQ >= 30
            if not _genotype_check(gt, gq):  
                continue
            
            altid = int(gt.split('/')[1])
            var = rec.chrom + ':' + str(rec.pos) + ':' + rec.ref +\
            ':' + rec.alt[altid - 1]
            af, flag = isRare(altid, rec.info, thres_af)
            if not flag:  # Checks if variant is not Rare (AF < 5%) in ExAC
                continue
            
            if 'LCR' in rec.info:
                continue
            
            if 'CLNDBN' in rec.info:
                dn = rec.info.CLNDBN[altid - 1]
                sig_num = rec.info.CLNSIG[altid - 1]
                if '|' in sig_num:
                    sig_num = [int(e) for e in sig_num.split('|') if e != '.']
                    if sig_num:
                        sig_num.sort()
                        sig_num = sig_num[-1]
                        cln_sig = CLNSIG_MAP[sig_num]
                    else:
                        cln_sig, dn = '', ''
                elif sig_num != '.':
                    sig_num = int(sig_num)
                    cln_sig = CLNSIG_MAP[sig_num]
                else:
                    dn, cln_sig = '', ''
            else:
                dn, cln_sig = '', ''
                
            if 'LCR' in rec.info:
                lcr = 'LCR'
            else:
                lcr = ''
            
            if 'CADD_phred' in rec.info:
                val = rec.info['CADD_phred'][altid - 1]
                if val == '.':
                    cadd = ''
                else:
                    cadd = float(val)
            else:
                cadd = ''
            
            if len(rec.ref) == len(rec.alt[altid - 1]) and len(rec.ref) == 1:
                ada_score, rf_score = get_dbscSNV_ant(rec.chrom, rec.pos,
                                                rec.ref, rec.alt[altid - 1])
                if (ada_score and ada_score > 0.6) or (rf_score and
                                                       rf_score > 0.6):
                    scpred = 'Damaging'
                else:
                    scpred = ''
            else:
                ada_score, rf_score, scpred = '', '', ''
            sc_ant = [scpred, ada_score, rf_score]
            
            #Parse annotation and prioritize transcript
            pa = vp.prio_trans(vp.parse(rec.info))
            
            # Ignore the intergenic variants
            if altid not in pa:  
                continue
            
            eqtl_flag = False
            for gene, ant in pa[altid].items():
                ta = ant['TRANSCRIPT']
                key = ta.trans_id + '_' + ta.aa
                snps3d_pred = ['', '', '', '']
                
                # SC-1 variant present in Clinvar as Pathogenic or Likely Pathogenic
                if sc == 1: # Search Criteria 1
                    if (cln_sig in ['Pathogenic', 'Likely pathogenic']):
                        if gene not in data[sid]:
                            data[sid][gene] = []
                        data[sid][gene].append((ta, af, cadd, gt, eqtl_flag,
                                                var, dn, cln_sig, lcr,
                                                sc_ant, snps3d_pred))
                                        
                # SC-2 variant is protein altering + SC-1
                if sc == 2: # Search Criteria 2
                    if (_is_PASnv(ta) and _is_Damaging(altid, rec.info, ta,
                                snps3d_pred, nmethod)) or _is_NonSense(ta) \
                                or _is_Splicing(ta) or _is_PAIndel(ta) or \
                                scpred == 'Damaging' or (cln_sig in ['Pathogenic',
                                               'Likely pathogenic']):
                        if gene not in data[sid]:
                            data[sid][gene] = []
                        data[sid][gene].append((ta, af, cadd, gt, eqtl_flag,
                                                var, dn, cln_sig, lcr,
                                                sc_ant, snps3d_pred))
                
                # SC-5 in intronic and UTR variants + SC-1 + SC-2
                if sc == 3: # Search Criteria 3 
                    if _is_Intronic(ta) or _is_UTR(ta) or _is_PASnv(ta) or \
                    _is_NonSense(ta) or _is_Splicing(ta) or _is_PAIndel(ta) \
                    or scpred == 'Damaging' or cln_sig in ['Pathogenic',
                    'Likely pathogenic']:
                        if gene not in data[sid]:
                            data[sid][gene] = []
                        data[sid][gene].append((ta, af, cadd, gt, eqtl_flag,
                                                var, dn, cln_sig, lcr,
                                                sc_ant, snps3d_pred))
    return data


class VarPrioritization:
    
    def  __init__(self, genelist_config, gender_config):
        self.genelist_file = genelist_config
        self._load_genelist()
        self.gender_file = gender_config
        self._load_geneder()
        
    
    def _load_genelist(self):
        self.genelist = {}
        s = open(self.genelist_file)
        for line in s:
            line = line.strip()
            if not line:
                continue
            if '#' in line:
                disease_class = line[1:]
                if disease_class not in self.genelist:
                    self.genelist[disease_class] = {}
                continue
            if line.startswith('--'):
                continue
            d = line.split('\t')
            gene, exon, mod = d
            self.genelist[disease_class][gene] = (exon, mod)
        s.close()
                
    def _load_geneder(self):
        self.gender = {}
        s = open(self.gender_file)
        for line in s:
            line = line.strip()
            if not line:
                continue
            if '#' in line:
                gen = line[1:]
                continue
            d = line.split('\t')[0]
            self.gender[d] = gen

    def filter_by_inheritence(self, data, sid):
        fdata = {}
        gen = self.gender[sid]
        mut_gl = set(data.keys())
        for dc, ele in self.genelist.items():
            dis_gl = set(ele.keys())
            if len(dis_gl - mut_gl) == len(dis_gl):
                continue
            mg = dis_gl & mut_gl
            mg = list(mg)
            for g in mg:
                exon, mod = ele[g]
                # AR one variant per gene
                if mod in ['AR', 'P-AR'] and len(data[g]) == 1: 
                    ta, af, cadd, gt, eqtl_flag, var, dn, cln_sig, lcr,\
                    sc_ant, snps3d_pred = data[g][0]
                    a1, a2 = gt.split('/')
                    if a1 == a2:
                        if dc not in fdata:
                            fdata[dc] = {}
                        if g not in fdata[dc]:
                            fdata[dc][g] = data[g]
                # AR more than one variant per gene
                elif ('AR' in mod or 'P-AR' in mod) and len(data[g]) >= 1: 
                    vflag = False
                    for vi in data[g]:
                        #print vi
                        ta, af, cadd, gt, eqtl_flag, var, dn, cln_sig, lcr, \
                        sc_ant, snps3d_pred = vi
                        a1, a2 = gt.split('/')
                        if a1 == a2:
                            vflag = True
                            gflag = True
                            if dc not in fdata:
                                fdata[dc] = {}
                            if g not in fdata[dc]:
                                fdata[dc][g] = [vi]
                            else:
                                fdata[dc][g].append(vi)
                    # For Compound Heterozygous cases
                    if vflag == False:
                        if dc not in fdata:
                            fdata[dc] = {}
                        if g not in fdata[dc]:
                            fdata[dc][g] = data[g]
                # AD one or more than one variant per gene
                elif 'AD' in mod:  
                    if dc not in fdata:
                        fdata[dc] = {}
                    if g not in fdata[dc]:
                        fdata[dc][g] = data[g]
                # X-Recessive solved gender wise for one or more than one
                #variant per gene
                elif mod == 'XR':  
                    if gen == 'male':
                        if dc not in fdata:
                            fdata[dc] = {}
                        if g not in fdata[dc]:
                            fdata[dc][g] = data[g]
                    elif gen == 'female':
                        for vi in data[g]:
                            ta, af, cadd, gt, eqtl_flag, var, dn, cln_sig, lcr, \
                            sc_ant, snps3d_pred = vi
                            a1, a2 = gt.split('/')
                            if a1 == a2:
                                flag = True
                                if dc not in fdata:
                                    fdata[dc] = {}
                                if g not in fdata[dc]:
                                    fdata[dc][g] = [vi]
                                else:
                                    fdata[dc][g].append(vi)
        return fdata
                
                            
def _load_ethnicity(ethnicity_config):
    ETH = {'East Asian': 'EAS',
           'African': 'AFR',
           'European': 'EUR'}
    eth_data = {}
    flag = False
    s = open(ethnicity_config)
    for line in s:
        line = line.strip()
        if not line:
            continue
        if flag == True:
            sid, e = line.split(':')
            eth_data[sid.strip()] = ETH[e.strip()]
        if line.startswith('Following'):
            flag = True
    return eth_data     
                        
                        
def get_norm_pos(var):
    chrom, pos, r, a = var.split(':')
    pos = int(pos)
    if len(r) == len(a):
        spos = pos
        epos = pos
        ref, alt = r, a
        t = 'snp'
    elif len(r) > len(a):
        ref = r[1:]
        alt = '-'
        spos = pos + 1
        epos = pos + len(ref)
        t = 'del'
    elif len(r) < len(a):
        ref = '-'
        alt = a[1:]
        spos = pos
        epos = pos + 1
        t = 'ins'
    return chrom, spos, epos, ref, alt, t


  

def main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config,
         prev_patients, nmethod, thres_af = 0.00, sc = 1):
    
    vp = VarPrioritization(genelist_config, gender_config)
    
    ldata = {}
    ldata = _load(snv_vcf, thres_af, nmethod, ldata, sc)
    ldata = _load(indel_vcf, thres_af, nmethod, ldata, sc)
    
    eth_data = _load_ethnicity(ethnicity_config)
    
    if thres_af == 0:
        id = 'SC-%d' % sc + '_MAF=0'
    else:
        id = 'SC-%d' % sc + '_MAF<=%f' % thres_af
    c = 0
    for sid, data in ldata.items():
        fdata = vp.filter_by_inheritence(data, sid)
        if fdata:
            if sid in prev_patients:
                continue
            prev_patients.append(sid)
            c += 1
            for k, y in fdata.items():
                disease_class = k
                for g, val in y.items():
                    for v in val:
                        ta, af, cadd, gt, eqtl_flag, var, dn, cln_sig, lcr, \
                        hgmd_ant, sc_ant, snps3d_pred = v
                        if ta.mutation:
                            mech = ta.mutation
                        elif ta.splice:
                            mech = ta.splice
                        elif ta.region:
                            mech = ta.region
                        else:
                            mech = ''
                        od = [id] + [str(e) for e in [c, vp.gender[sid],
                            eth_data[sid], sid, var, gt, af, g,
                            vp.genelist[k][g][1], ta.trans_id,
                            ta.exon, ta.utr_signal, ta.cdna, mech, ta.aa,
                            ta.sift, ta.pp2, cadd]+ snps3d_pred +sc_ant + [
                            dn, cln_sig]] + hgmd_ant + [disease_class]
                        print '\t'.join(od)
            print '\n'
    #print thres_af, 'Pred_Method=', nmethod, 'PATIENT CLASSIFIED=', c
    return prev_patients
                


if __name__ == '__main__':
    basedir = '/moulthome/kunduk/CAGI2015/HopkinChallenge'
    snv_vcf = basedir + '/SNVs/4-Hopkins_clinical_panel_SNV_varant_refgene.vcf'
    indel_vcf = basedir + '/INDELs/4-Hopkins_clinical_panel_InDel_varant_refgene.vcf'
    genelist_config = basedir + '/genelist_config_vpaper.tsv'
    gender_config = basedir + '/gender_config.tsv'
    ethnicity_config = basedir + '/4-Hopkins_clinical_panel_SNV_predicted_ethnicity_KK.tsv'
    
    
    # Logic that was used in CAGI2015
    # C-1 Present in HGMD or Clinvar and annotated as DM, DP or Pathogenic
    prev_patients = []
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.00, sc = 1)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.005, sc = 1)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.01, sc = 1)
    
    # C-2 Predicted Damaging NonSyn, FrameShift/Non-FrameShift, Direct Splicing, Damaging by dbscSNVs, NonSense + C-1
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.00, sc = 2)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.005, sc = 2)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.01, sc = 2)
    
    # C-3 Intronic + UTR + C-2 + C-1
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.00, sc = 3)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.005, sc = 3)
    prev_patients = main(snv_vcf, indel_vcf, genelist_config, gender_config, ethnicity_config, prev_patients, nmethod=1, thres_af = 0.01, sc = 3)
    #print prev_patients
    
    
    
    
    
        
    


