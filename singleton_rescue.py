#!/usr/bin/env python

###############################################################
#
#                      Singleton Rescue
#
# Author: Nina Wang
# Date Created: July 5, 2016
###############################################################
# Function: 
# Written for Python 3.5.1
# 
# uid_dict(bamfile): create bam_read dictionaries and counts number of reads
# read_mode(field, bam_reads): return most common occurrence in a specified field of read (e.g. cigar, etc)
# create_aligned_segment(bam_reads, sscs, sscs_qual): Create consensus bam read
#
# Inputs: 
# 1. A position-sorted paired-end BAM file containing reads with a duplex tag 
#    in the header.  
#
# Description:
# Helper functions for single strand and duplex consensus making
#
###############################################################

import pysam # Need to install
import collections
import re
import array
from random import randint
from argparse import ArgumentParser

from consensus_helper import *
#from SSCS_maker import consensus_maker
from SSCS_median_fastq import consensus_maker, mismatch_pos


###############################
##         Functions         ##
###############################

def duplex_consensus(seq1, seq2, readlength):
    '''(str, str, int) -> str
    Return consensus of 2 sequences with N for variant bases.
    '''
    consensus = ''
    if seq1 == seq2:
        consensus = seq1
    else:
        for i in range(readlength):
            try:
                if seq1[i] == seq2[i]:
                    consensus += seq1[i]
                else:
                    consensus += 'N'
            except:
                consensus += 'N'
    
    return consensus


def rescue_duplex(read_tag, duplex_tag, singleton_dict, sscs_dict = None):
    '''(str, str, dict, dict) -> Pysam.AlignedSegment
    
    Return consensus read from singleton + duplex read (either found in SSCS or singleton bam).
    - Duplex read: Matching sequence on the other strand
    
    Quality scores: Retain quality score of singleton being rescued
    
    If additional dict is provided (SSCS_dict), duplex rescue with SSCS_dict. If no additional bamfile provided, singleton_dict will bee used for read rescue.
    '''
    
    read = singleton_dict[read_tag][0]
    
    # If bamfile provided, perform SSCS singleton rescue 
    if sscs_dict == None:
        duplex_read = singleton_dict[duplex_tag][0]
    else:
        duplex_read = sscs_dict[duplex_tag][0]
        
    dcs = duplex_consensus(read.query_alignment_sequence, duplex_read.query_alignment_sequence, read.query_alignment_length)

    # If consensus has >30% N's, toss read
    if dcs.count('N')/len(dcs) > 0.3:
        return None
    
    dsc_read = create_aligned_segment([read], dcs, read.query_alignment_qualities)
    
    ## If SSCS used to rescue read, create segment based on 
    #if bamfile == None:
        #dsc_read = create_aligned_segment([read], dcs, read.query_alignment_qualities)
    #else:
        #dsc_read = create_aligned_segment([read, duplex_read], dcs, read.query_alignment_qualities)
        

    return dsc_read


def pileup(pysam_bam, read_chr, start, read_cutoff):
    '''(pysam.calignmentfile.AlignmentFile, str, int, int) -> str
    
    Return consensus of pileup reads at a single position. 
    
    - Read_cutoff: 1 SSCS or 3 Singleton (as read is included in pileup) => 1.0 allelic freq cutoff (want it to be more stringent since reads are not from same mol)
    - Should we make exceptions for this? e.g. [('T', 4812), ('N', 1)] => N due to soft clips
    
    Return None or 'N'???
    '''    
    for pileupcol in pysam_bam.pileup(read_chr, start, start + 1, truncate = True):
        #nuc = [0, 0 ,0, 0, 0] # A, C, G, T, N 
        nuc = ''
        
        # Make consensus from > 2 reads at each position
        if pileupcol.n < read_cutoff:
            return 'N'
        else:
            #print(pileupcol)
            for pileupread in pileupcol.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:                       
                    nuc += pileupread.alignment.query_sequence[pileupread.query_position]
                    #if pileupread.alignment.query_sequence[pileupread.query_position] == 'N':
                        #print(pileupread)
                        #print(pileupread.query_position)

                    
        consensus = collections.Counter(nuc).most_common(1)[0]  
        
        # Since we trust pileup reads less, set more stringent threshold (all the reads have to have same base to make consensus)
        #print(consensus)
        #print(len(nuc))
        #print(collections.Counter(nuc).most_common())
        if consensus[1]/len(nuc) != 1:
            return 'N'
        else:
            return consensus[0]    


def rescue_pileup(pysam_bam, bam_read, read_cutoff):
    '''(str, pysam.calignmentfile.AlignmentFile, pysam.calignedsegment.AlignedSegment, int) -> str
    
    Return consensus sequence from pileup of all reads in region of bam_read (singleton).
    
    DETERMINE LOCATION OF VARIANT -> if less than specified # of reads, don't make consensus (utilize SSCS_maker.py - mismatch_pos)
    
    Require variant positions to meet read cutoff rather than each position of sequence since non-variant positions match to reference and shouldn't differ by much. (Also faster runtime if we only do mpileup at variant positions!)
    - Read_cutoff: 1 SSCS or 3 Singleton (as read is included in pileup) => 1.0 cutoff (want it to be more stringent since reads are not from same mol)
    
    If SSCS pileup rescue, compare SSCS consensus with singleton being rescued => mismatch positions assign N
    '''
    start = bam_read.reference_start
    
    consensus = ''
    
    # === Inorporate soft clipped regions with start and stop pos ====
    readLength = bam_read.query_alignment_length # soft clips not included
    
    if 'S' in bam_read.cigarstring:
        cigar_pos = [x for x,y in bam_read.cigar] # => IT COULD BE BOTH FRONT AN END
        soft_pos = [i for i,x in enumerate(cigar_pos) if x == 4]
        
        if len(soft_pos) == 1:
            # Number of soft clips at start
            start_sc = bam_read.cigar[soft_pos[0]][1]
                
        else:
            start_sc = bam_read.cigar[soft_pos[0]][1]
            # Number of soft clips at end
            end_sc = bam_read.cigar[soft_pos[1]][1]
                        
        consensus += 'N'*start_sc
            
    
    mismatch_pos_lst = mismatch_pos(bam_read.cigar, bam_read.get_tag('MD'))
    
    for i in range(readLength):
        if i in mismatch_pos_lst:
            # Pileup consensus of variant
            mismatch_start = start + i
            
            consensus += pileup(pysam_bam, bam_read.reference_name, mismatch_start, read_cutoff)
            
        else:
            # No pileup needed for bases matching ref
            consensus += bam_read.query_alignment_sequence[i] # soft clips not included
            
    if 'end_sc' in locals():
        consensus += 'N'*end_sc

    if read_cutoff == 1:
        # Compare singleton read with SSCS consensus (not needed for singleton consensus as singleton already included)
        consensus = duplex_consensus(bam_read.query_sequence, consensus, bam_read.query_length)
    
    consensus_read = create_aligned_segment([bam_read], consensus, bam_read.query_alignment_qualities)
    # Check to see consensus < 0.3 Ns (excluding softclips)
    if consensus.count('N')/len(dcs) > 0.3:
        return None
    
    return consensus
    
    
def main():
    # Command-line parameters
    parser = ArgumentParser()
    parser.add_argument("--infile", action = "store", dest="infile", help="input singleton BAM file", required = True)
    #parser.add_argument("--rescue_outfile", action = "store", dest="rescue_outfile", help="output BAM file", required = True)
    args = parser.parse_args()
    
    #infile ='/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test/MEM-001_KRAS.singleton.sort.bam'

    ### Later include flags to turn on/off certain rescues ###
    
    
    start_time = time.time()
    
    # initialize bams and dicts
    singleton_bam = pysam.AlignmentFile(args.infile, "rb") 
    sscs_bam = pysam.AlignmentFile('{}.sscs.sort.bam'.format(args.infile.split('.singleton.sort.bam')[0]), "rb")
    rescue_bam = pysam.AlignmentFile('{}.rescue.bam'.format(args.infile.split('.singleton.sort.bam')[0]), "wb", template = sscs_bam)
    stats = open('{}_stats.txt'.format(args.infile.split('.singleton.sort.bam')[0]), 'a')
    time_tracker = open('{}_time_tracker.txt'.format(args.infile.split('.singleton.sort.bam')[0]), 'a')
    
    singleton_dict = uid_dict(singleton_bam)[0]
    sscs_dict = uid_dict(sscs_bam)[0]
    
    
    sscs_dup_rescue = 0
    singleton_dup_rescue = 0
    sscs_rescue = 0
    singleton_rescue = 0
    
    rescue_dict = collections.OrderedDict()
    
    for i in singleton_dict.keys():
        barcode = i.split('_')[0]
        barcode_bases = int(len(barcode)/2)
        pair_barcode = barcode[barcode_bases:] + barcode[:barcode_bases] # get barcode of mate pair read (opposite strand)
        
        read_num = i[-2:]
        if read_num == 'R1':
            read_num = 'R2'
        else:
            read_num = 'R1'
        
        duplex = pair_barcode + '_' + i.split('_', 1)[1][:-2] + read_num # pair read is other read on opposite strand (e.g. read = ACGT_fwd_R1, duplex = GTAC_fwd_R2)
        
        
        # Check if there's duplicate singletons or if singleton is already rescued (in the case of singleton - singleton duplex rescue) => not an issue if we're not rescuing singleton & duplex at same time
        #if i not in rescue_dict.keys():
        
        # 1) Duplex rescue: Singleton + SSCS -> Check for singleton matching duplex sequence in SSCS bam
        if duplex in sscs_dict.keys():
            rescue_read = rescue_duplex(i, duplex, singleton_dict, sscs_dict)
            # If read doesn't passes Ncutoff (0.3), None is given
            if rescue_read != None:
                sscs_dup_rescue += 1
                
                #dcs_read.query_name = "{}|{}|{}\t".format(dcs_read.query_name.split('|')[0], min(barcode, pair_barcode), max(barcode, pair_barcode)) # Add both barcodes to header in alphabetical order [NOT CURRENTLY IMPLEMENTED - as rescued reads still go through duplex consensus making]
        
        # 2) Duplex rescue: Singleton + Singleton -> Check for singleton matching duplex sequence in Singleton bam -> WOULD YOU RESCUE BOTH SINGLETONS THEN??
        elif duplex in singleton_dict.keys():
            rescue_read = rescue_duplex(i, duplex, singleton_dict)
            
            if rescue_read != None:
                singleton_dup_rescue += 1                
        
        # 3) Non-duplex rescue: Any SSCS mapping to region => use pileup method            
        else:
            pileup_consensus = rescue_pileup(sscs_bam, singleton_dict[i][0], 1)
            
            if pileup_consensus != None:
                sscs_rescue += 1
            else:
                # 4) Non-duplex rescue: Any singleton mapping to region => use pileup method 
                pileup_consensus = rescue_pileup(singleton_bam, singleton_dict[i][0], 3)
                
                if pileup_consensus != None:
                    singleton_rescue += 1
                    
            
                
                   
                ##### Hold off on singleton_rescue right now as singleton seq soft clips are diff -> soft clips not included in sequence (so shorter seq, need to account for that) 
                
                #readLength = max(collections.Counter(i.query_alignment_length for i in singleton_dict[i]))
                #sscs_consensus_seq = rescue_pileup(i, singleton_bam, readLength, 3)
                #if sscs_consensus_seq != None:
                    #singleton_rescue += 1                    

            try:
                rescue_read = create_aligned_segment([singleton_dict[i][0]], sscs_consensus_seq, singleton_dict[i][0].query_alignment_qualities)
            except:
                #print(readLength)
                print(sscs_consensus_seq)
                print(singleton_dict[i][0])             
                return 'hi'
                # Note: rescued seq does not include soft clips 
                # e.g. GCCT_chr12_25398127_25398127_rev_R2
                # singleton_seq = NNNNNNNNNNNNNATGAAAATGGTCAGAGAAACCTTTATCTGTATCAAAGAATGGTCCTGCACCAGTAATATGCATATTAAAACAAGATTTACCTCTA
                # consensus = ATGAAAATGGTCAGAGAAACCTTTATCTGTATCAAAGAATGGTCCTGCACCAGTAATATGCATATTAAAACAAGATTTACCTCTATTGTTGGATCATA
                
                #print(sscs_consensus_seq)
                #if sscs_consensus_seq != None:
                    #sscs_rescue += 1
                    #read = singleton_dict[i][0]
                    #rescue_read = create_aligned_segment([read], sscs_consensus_seq, read.query_alignment_qualities)
                    
                #else:
                    #singleton_consensus_seq = rescue_pileup(i, singleton_bam)
                    #if singleton_consensus_seq !=None:
                        #singleton_rescue += 1
                        #read = singleton_dict[i][0]
                        #rescue_read = create_aligned_segment([read], singleton_consensus_seq, read.query_alignment_qualities)               
            
            
            if rescue_read != None:
                rescue_dict[i] = singleton_dict[i][0]
                rescue_bam.write(rescue_read)
            
    time_tracker.write('Singleton Rescue: ')
    time_tracker.write(str((time.time() - start_time)/60) + '\n')

    
    summary_stats='''SSCS duplex rescue: {} \n
Singleton duplex rescue: {} \n
SSCS rescue: {} \n
Singleton rescue: {} \n'''.format(sscs_dup_rescue, singleton_dup_rescue, sscs_rescue, singleton_rescue)
    
    stats.write(summary_stats)
    
    
    singleton_bam.close()
    sscs_bam.close()
    rescue_bam.close()
    stats.close()
    time_tracker.close()
    

###############################
##           Main            ##
###############################
if __name__ == "__main__": 
    import time
    start_time = time.time()
    main()  
    
    ## TEST
    #infile = '/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/pl_duplex_sequencing/test2/MEM-001_KRAS.singleton.sort.bam'
    #singleton_bam = pysam.AlignmentFile(infile, "rb")
    #sscs_bam = pysam.AlignmentFile('{}.sscs.sort.bam'.format(infile.split('.singleton.sort.bam')[0]), "rb")

    #bam = pysam.AlignmentFile('/Users/nina/Desktop/Ninaa/PughLab/Molecular_barcoding/code/PL_consensus_test/test2/MEM-001_KRAS.bam', "rb")
    #bam_dict = uid_dict(bam)[0]
    
    #singleton_dict = uid_dict(singleton_bam)[0]
    #sscs_dict = uid_dict(sscs_bam)[0]    
    
    ## Find key containing barcode and start coor
    # [x for x in singleton_dict.keys() if re.match('ATCG_chr12_25398169', x)]
    
    #i = 'GTGG_chr12_25398271_25398214_rev_R1'
    #barcode = i.split('_')[0]
    #barcode_bases = int(len(barcode)/2)
    #pair_barcode = barcode[barcode_bases:] + barcode[:barcode_bases] # get barcode of mate pair read (opposite strand)
    #read_num = i[-2:]
    #if read_num == 'R1':
        #read_num = 'R2'
    #else:
        #read_num = 'R1'
    
    #duplex = pair_barcode + '_' + i.split('_', 1)[1][:-2] + read_num # pair read is other read on opposite strand (e.g.
    
    ##rescue_pileup(i, sscs_bam)
    #readLength = sum([y for x,y in singleton_dict[i][0].cigar])
    
    #rescue_pileup('GCCT_chr12_25398127_25398127_rev_R2', sscs_bam, readLength)


    print((time.time() - start_time)/60)  
            
            
            
            #elif duplex in singleton_dict.keys() and duplex not in rescue_dict.keys():
                #rescue_read = rescue_duplex(i, duplex, singleton_dict)
                
                #if dcs_read != None:
                    #singleton_dup_rescue += 1
                    ##dcs_read.query_name = "{}|{}|{}\t".format(dcs_read.query_name.split('|')[0], min(barcode, pair_barcode), max(barcode, pair_barcode)) # Add both barcodes to header in alphabetical order
                    #rescue_bam.write(rescue_read)
                    
                    ## rescue duplex pair
                    #duplex_rescue_read = rescue_duplex(duplex, i, singleton_dict)
                    #rescue_bam.write(duplex_rescue_read)
            
            ## 3) Non-duplex rescue: Any SSCS mapping to region => use pileup method
            #else:
                #sscs_consensus_seq = rescue_pileup(i, sscs_bam)
                
                #if sscs_consensus_seq != None:
                    #sscs_rescue += 1
                    #read = singleton_dict[i][0]
                    #rescue_read = create_aligned_segment([read], sscs_consensus_seq, read.query_alignment_qualities)
                    
                #else:
                    #singleton_consensus_seq = rescue_pileup(i, singleton_bam)
                    #if singleton_consensus_seq !=None:
                        #singleton_rescue += 1
                        #read = singleton_dict[i][0]
                        #rescue_read = create_aligned_segment([read], singleton_consensus_seq, read.query_alignment_qualities)                        
            
            #if rescue_read != None:
                #rescue_bam.write(rescue_read)
            
            

            
            
            #read_chr = i.split('_')[1]
            #read_start = i.split('_')[2]
            #read_end = i.split('_')[3]            
            
            
            #read_fwd = '{}_{}_{}'.format(read_chr, read_start, read_end)
            #read_rev = '{}_{}_{}'.format(read_chr, read_end, read_start)

            #keys_in_region = [key for key in sscs_dict if read_fwd in key or read_rev in key]
            
            #if len(keys_in_region) > 0:
                #reads_in_region = [sscs_dict.get(key)[0] for key in keys_in_region]
                #reads_in_region += [singleton_dict.get(i)]
                
                #readLength = max(collections.Counter(read.query_alignment_length for read in reads_in_region))
                
                ## consensus maker to find most freq between SSCSs + singleton read
                #rescue_consensus = consensus_maker(reads_in_region, readLength, 0.7)
                
                #if rescue_consensus[0].count('N')/len(rescue_consensus[0]) > float(0.3):
                    #continue
                #else:
                    #rescue_read = create_aligned_segment([reads_in_region], rescue_consensus[0], rescue_consensus[1])
                    #rescue_bam.write(rescue_read)

                
            ## 4) Any read in region of variant
            #bam_in_region = sscs_bam.fetch(read_chr, int(read_start), int(read_end))
            #reads_in_region = [read for read in bam_in_region]
            
            
            #elif any('{}_{}_{}'.format(read_chr, read_start, read_end) in s for s in sscs_dict.keys()) or any('{}_{}_{}'.format(read_chr, read_end, read_start) in s for s in sscs_dict.keys()):
                
                
                
            
            

            #for read in sscs_bam.fetch(read_chr, int(read_start), int(read_end)):
                
            
            
            
            #else:

# rescued bams should be based on the singleton sequence, but with the duplex consensus 
# duplex consensus should have filters (Ncutoff) to reflect 



#### pileup_rescue 

#start = read_start
#readLength = sum([y for x,y in bam_read.cigar])
#end = int(read_start) + readLength

#### NEED TO ADDRESS SOFT CLIPS AT START OF SEQ -> actual seq much shorter
## don't locate with start and end, use start and read length to determine consensus 
##if (strand == 'fwd' and read == 'R1') or (strand == 'rev' and read == 'R1'):
    ##start = read_start
##else:
    ##start = read_end

##end = int(start) + readLength

 
##else:
    ##if strand == 'fwd':
        ##end = read_end
    ##else:
        ##end = read_start
    ##start = int(end) - readLength        

#consensus = ''

### Soft clips: only want seq that match ref, give Ns for remainder => mpileup shorter region to reflect what is sequenced

##sc_end = False
###readLength = sum([y for x,y in read.cigar])
#if 'S' in bam_read.cigarstring:
    #soft_pos = [x for x,y in bam_read.cigar].index(4) # => IT COULD BE BOTH FRONT AN END
    #num_sc = bam_read.cigar[soft_pos][1]
    
    #end = int(read_start) - readLength
    ##if soft_pos == 0:
        ##consensus += 'N'*num_sc
    ##else:
        ##sc_end = True
        ##start = int(read_start) + num_sc
    ##else:
        ##end = int(read_start) + readLength    

##print(start,end)
#### PILE UP ONLY VARIANT REGIONS 
#for pileupcol in pysam_bam.pileup(read_chr, int(start), int(end)):
    #nuc = '' # change to dictionary
    ##print(pileupcol.n)
    ## Make consensus from > 2 reads at each position
    #if pileupcol.pos < int(start) or pileupcol.pos >= int(end):
        #continue
    #elif pileupcol.n < read_cutoff: ## <== SHOULD WE RESTRICT THIS TO VARIANT POS ONLY?? yes! determine variant from singleton bam
        #return None
    #else:
        ##print(pileupcol)
        #for pileupread in pileupcol.pileups:
            #if not pileupread.is_del and not pileupread.is_refskip:

                #try:
                    ##print(pileupread.query_position)
                    ##print(pileupread.alignment.query_sequence)
                    ##print(pileupread.alignment.query_sequence[pileupread.query_position])  
                    ##print(pileupread)                        
                    #nuc += pileupread.alignment.query_sequence[pileupread.query_position]
                #except:
                    #print(pileupread.query_position)
                    #print(pileupread.alignment.query_sequence)
                    #print(pileupread.alignment.query_sequence[pileupread.query_position - 1])  
                    #print(pileupread)
                    ##print(pileupread.next_reference_start)
                    ##print(pileupread.is_reverse)
                    ##print(pileupread.is_read2)
                    
                    
        ##print(collections.Counter(nuc).most_common(1)[0][0])
        ##return collections.Counter(nuc).most_common(1)
        
        #### INCLUDE SINGLETON WITH SSCS reads => after making consensus, compare consensus with singleton (0.7 cut off for variants) 
        #consensus += collections.Counter(nuc).most_common(1)[0][0]

##if sc_end == True:
    ##consensus += 'N'*num_sc
#return consensus

