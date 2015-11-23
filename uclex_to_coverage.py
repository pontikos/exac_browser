import sys
import os
import pysam
from pysam import VCF
import numpy
import csv
import cPickle as pickle


BP_file =  "BP_v_GEN19_COMP_sDaiu_analysis_ready.sort.bed" #sys.argv[1]

#gvcf_chromosome = "/SAN/biomed/biomed5/biomed5/ExAC_data/version2/ExAC.r0.2.sites.vep.vcf.gz" #sys.argv[1] #"/SAN/biomed/biomed6/biomed6/warren/BIG_vcf/Chr22_02.recode.vcf.gz" #"/SAN/biomed/biomed6/biomed6/warren/BIG_vcf/chr22_005.vcf.gz" #"/SAN/biomed/biomed6/biomed6/warren/BIG_vcf/Chr22_02.recode.vcf.gz" #"/SAN/biomed/biomed6/biomed6/warren/BIG_vcf/Chr22_02.recode.vcf.gz" #"/SAN/biomed/biomed6/biomed6/warren/BIG_vcf/chr22_test2.recode.vcf.gz" #  "chr22_test2.recode.vcf.gz" # "chr22_test2.recode.vcf.gz" #sys.argv[1] # "chr22_test2.recode.vcf.gz" #
vcf_file =  sys.argv[1]


#type_ = os.path.split(gvcf_chromosome)[1]

chr_ = vcf_file.split(".")[1]


name = "%s_%s"%(BP_file.split(".")[0],chr_)

'''
# take branchpoint file

1. calculate a area around it, use this to plot the graph
2. find variants overlapping position 1/3 ----- also composition of transition/transversions and changes A>N A>indel
3. split data into subcatagories and output relevant files

'''

    
############## Parameters
flanking_region = 20
total_length = flanking_region*2+1


#positions of interest, 0 being the middle position, -2 being 2 nucleotides upsteam : POSTIVE IMPLIES UPSTREAM
positions_of_interest_index = [0,2]
positions_of_interest_names = ['Pos_3', 'Pos_1']
#"""

position_of_exon_bound = 3
   

######## Csv input

ifile_bp  = open(BP_file, "rb")
readerbp = csv.reader(ifile_bp,delimiter='\t')

############ Dictionaries

## capture variation at key locations
variant_dict = {}


## graph dictionaries

# overall 
graph_dict = {"ALL":{}}

for group_type in graph_dict:
    
    for i in range(1,total_length+1):

        graph_dict[group_type][i] = {"variant":0,"samples":0} 


# overall multi-count
graph_dict_full_counts = {"ALL":{}}

for group_type in graph_dict_full_counts:
    
    for i in range(1,total_length+1):

        graph_dict_full_counts[group_type][i] = {"variant":0,"samples":0} 
            
            
    

#############
### WARNING - VCF encoded at 1-base and BED encoded at 0 base, pysam converts to 0-based...
#############
#'''
vcf = VCF()

vcf.connect(vcf_file)

samples = vcf.getsamples()
#'''


ofile  = open('%s_variants.csv'%(name), "wb") # ofile  = open("chr22_c2_pass14_class_BP.csv", "wb")
writer_variants = csv.writer(ofile, quoting=csv.QUOTE_NONE)


#combined_BP_dict = {}
count = 0
for bp_region in readerbp:
    
    
    
    #raw_input()

    [chr_bp,bp_left,bp_right,bp_id,x,strand_bp] = bp_region[:6]
    
    if chr_bp != chr_:
        continue
    
    #print bp_region
    
    #print chr_bp
    #print chr_
    #count +=1
    #raw_input()
    BP_type = "ALL"
    
    
    # splitting bp_id 
    #  'NEJC_chr1_985588_985589|ENST00000379370.2_AGRN|1'
    #  'MERC_chr1_981087_981088_C|ENST00000379370.2_AGRN|4'
    #'  BOTH_ENSG00000188157__chr1_981081_981082_A|ENST00000379370.2_AGRN|4'
    #print bp_id
    #[ens_id,gene_name]=bp_id.split("|")[1].split("_")
    #[bp_intron_count]=bp_id.split("|")[2]
       
   
    #print [chr_bp,bp_left,bp_right,bp_id,x,strand_bp] 
    
    chr_bp = chr_bp.replace("chr","")
    
    bp_left = int(bp_left)
    bp_right = int(bp_right)
    
    positions_of_interest = []

    #print poi
    #print bp_right
    
    #'''
    if strand_bp == "+":

        for poi_i in positions_of_interest_index:
                positions_of_interest.append(bp_right- poi_i)            
          
    if strand_bp == "-":
        for poi_i in positions_of_interest_index:
            positions_of_interest.append(bp_right + poi_i)
    
    
    #'''        
    
    
    # create window around BP 3rd pos
    
    left_window = bp_right - flanking_region
    right_window= bp_right+ flanking_region
    
    #print "\t".join(map(str,[chr_bp,bp_left,bp_right,bp_id]))
    #print "\t".join(map(str,[chr_bp,left_window,right_window,"Windw"]))
    
    #print "left_window",left_window
    #print "right_window",right_window      
    
    #raw_input()
    #select out start / end dependant on the strand and position of exon!!!!!!
    bp_start = 0
    #'''
    if strand_bp == "+":        
        bp_start = left_window
    else: 
        bp_start = right_window
    #'''
    #bp_start = bp_right
    
    #print "bp_start",bp_start
    #print "abs(bp_start-positions_of_interest[0])",abs(bp_start-positions_of_interest[0])        
    #raw_input()
    
    #vcf_bp = vcf.fetch(chr_bp,left_window,right_window)    
    try:
        #vcf_bp = vcf.fetch("22",17619278,17619309)
        vcf_bp = vcf.fetch(chr_bp,left_window,right_window)
    except:
            
        print "Failed to load region: ",chr_bp,left_window,right_window
        continue    
    
    pos_in_list = range(1,total_length+1)
    
    found_pos = {}
    
    
    
    # for sample,record,geno in zip(samples,vcf_bp,record_data):
    
    
    # get all variants within the region around splice site
    for record in vcf_bp:
        #for record in vcf_bp:
        #print record.info
        #print record.filter
        #print record.qual
        
                
        # save homozygous samples 
        sample_homs = []
        
        # save hets samples 
        sample_hets = []        
        
        ## counts
        count_dict_homo = 0  
        count_dict_hetro = 0
        total_calls = 0            
        
        position = record.pos +1
        
        
        
        if position > left_window and position < right_window:
            
            # relative to window around branchpoint
            relative_pos=  abs(bp_start-position)
            
            record_data = str(record).split("\t")[9:]
            
            # remove position 
            pos_in_list.remove(relative_pos)

            #print position
        
            # single count add sample
            graph_dict["ALL"][relative_pos]["samples"] += 1   
            
            
            
            ################  do all the counts
            for sample,geno in zip(samples,record_data):
            
                gt = geno.replace("|","/").split(":")[0].split("/")
                
                if gt == ["."]:
                                       
                    #print "NO GT, skipping..."
                    #raw_input()
                    continue        
        
                # count
                total_calls+=1
                graph_dict_full_counts["ALL"][relative_pos]["samples"] +=1
                
                if gt.count("0") == 0:
                    #print "homo!!! :",record[sample]['GT'][0]
                    count_dict_homo +=1
                    #graph_dict[BP_type][relative_pos]["variant"] += 1
                    graph_dict_full_counts["ALL"][relative_pos]["variant"] += 1            
                    #raw_input("FULL HOMO")
                    sample_homs.append([sample,"|".join(gt)])
                
                elif gt.count("0") == 1:
                    count_dict_hetro +=1      
                    #graph_dict[BP_type][relative_pos]["variant"] += 1
                    graph_dict_full_counts["ALL"][relative_pos]["variant"] += 1        
                    
                    sample_hets.append([sample,"|".join(gt)])


            # adding 1 to single count data
            if count_dict_hetro+count_dict_homo >= 1:

                graph_dict["ALL"][relative_pos]["variant"] += 1
                
            
            
            # if variant falls on positions of interest ie. possible deleterious variant
            if position in positions_of_interest:
    
                relative_bp_loc = positions_of_interest_names[positions_of_interest.index(position)]
                #print "alignment check"
                
                #print strand_bp
                
                
                '''
                print " ".join(bp_region)
                print record.ref
                print record.alt
                #print record
                
                print position
                
                print relative_bp_loc
                
                #print positions_of_interest
                             
                #print positions_of_interest[0]
                #print positions_of_interest[position_of_exon_bound]
                
                print "bp_left",bp_left
                print "bp_right",bp_right
                raw_input()  
                #'''
                
                for samp_gt in sample_homs:
                    writer_variants.writerow(samp_gt+["HOM",record.contig,record.id,position,bp_id,strand_bp,relative_bp_loc,record.ref,record.alt,bp_right,total_calls,count_dict_hetro,count_dict_homo])
            
                for samp_gt in sample_hets:
                    writer_variants.writerow(samp_gt+["HET",record.contig,record.id,position,bp_id,strand_bp,relative_bp_loc,record.ref,record.alt,bp_right,total_calls,count_dict_hetro,count_dict_homo])                
                
    # add missed positions with no variants!!
    if len(pos_in_list) > 0:
        
        for missed_pos in pos_in_list:
    
            # single counts
            graph_dict["ALL"][missed_pos]["samples"] += 1   
            #full counts
            #graph_dict_full_counts[BP_type][missed_pos]["samples"] = len(samples)
            graph_dict_full_counts["ALL"][missed_pos]["samples"] += len(samples)              
            
    
 
        
       
              

pickle.dump(graph_dict,open("%s_graph_dict.p"%(name),"w"))
pickle.dump(graph_dict_full_counts,open("%s_graph_dict_full_counts.p"%(name),"w"))
#"""

#graph_dict = pickle.load(open("/Users/warren/scripts/Branchpoints/Nejc_data/gVCF/fullhom_jan2015/EXAC/results/BP_v_GEN19_COMP_sDaiu_analysis_ready_graph_dict.p"))
###### print graphs

# graph types

#graphs_list = [graph_dict,graph_dict_full_counts]


ofile4  = open('%s_COUNT_BP.csv'%(name), "wb") # ofile  = open("chr22_c2_pass14_class_BP.csv", "wb")
writer_count = csv.writer(ofile4, quoting=csv.QUOTE_NONE)    

ofile3  = open('%s_RATIO_BP.csv'%(name), "wb") # ofile  = open("chr22_c2_pass14_class_BP.csv", "wb")
writer_ratio = csv.writer(ofile3, quoting=csv.QUOTE_NONE)

files = []
for graph_type in graph_dict:

    
    
  

    for pos_ in range(1,total_length+1):
        
        if float(graph_dict[graph_type][pos_]["samples"]) == 0.0:
            
            print "no count for pos_!!"
        
        else:
            writer_ratio.writerow(["Splice_site_5p",graph_type,pos_,float(graph_dict[graph_type][pos_]["variant"])/float(graph_dict[graph_type][pos_]["samples"])])
            
        
        writer_count.writerow(["Splice_site_5p",graph_type,pos_,float(graph_dict[graph_type][pos_]["variant"])])

        




print count

