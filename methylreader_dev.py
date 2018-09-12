#this is to reade methylation patterns from a psecific targeted region##

# Author Arna

import sys
import pysam
import os
import re

from argparse import ArgumentParser



parser = ArgumentParser()

parser.add_argument("-b", "--bamfile", dest="bam",  help="The bam file output from bismark", type=str)

parser.add_argument("-L", "--interval", dest="interval",  help= "The genomic interval", type=str)

parser.add_argument("-wd", "--working_dir", dest="cwd",  help="Working directory where all the bam files are ", type=str)

parser.add_argument("-o", "--output_dir",  dest="output", help="Output directory", type=str)

args = parser.parse_args()

file = args.bam
target = args.interval
cdir = args.cwd
outputdir = args.output

print "Input bam is", file
print "Interval is ", target
print "Output direcotry is ", outputdir
print "Working directory is ", cdir



#######################					  
#     functions 	
#######################


def top_n(word_list, n):
    word_counter = {}
    for word in word_list:
        if word in word_counter:
         word_counter[word] += 1
        else:
            word_counter[word] = 1

    popular_words = sorted(word_counter, key = word_counter.get, reverse = True)
 
    top_nn = popular_words[:n]
    return top_nn

def unmethylpos(string):

	unmethylPos= [x.start() for x in re.finditer('z', string)]
	return unmethylPos	

def methylpos(string):

	methylPos =  [x.start() for x in re.finditer('Z', string)]

	return methylPos

def unmethylCitocineCHH(string):

        unmethylCitocineCHH =  [x.start() for x in re.finditer('h', string)]

        return unmethylCitocineCHH

def methylCitocineCHH(string):

        methylCitocineCHH =  [x.start() for x in re.finditer('H', string)]

        return methylCitocineCHH


def unmethylCitocineCHG(string):

        unmethylCitocineCHG =  [x.start() for x in re.finditer('x', string)]

        return unmethylCitocineCHG

def methylCitocineCHG(string):

        methylCitocineCHG =  [x.start() for x in re.finditer('X', string)]

        return methylCitocineCHG




def addLocus(sites,pos):
        newsites=[]
        for i in range(len(sites)):
                newsites.append(sites[i]+pos+1)
        return newsites
        
 
def subLocus(sites,pos):
        newsites=[]
        for i in range(len(sites)):
                newsites.append(pos - sites[i])
        return newsites  
        
             

def most_common(l):
    max = 0
    maxitem = None
    for x in set(l):
        count =  l.count(x)
        if count > max:
            max = count
            maxitem = x
    return maxitem
    

########

dir= os.path.dirname(file) 
string = re.split('[: -]',target)
print string
chrom = string[0]
bedS = int(string[1])
bedE = int(string[2])
bam = os.path.basename(file)


allReadsSites =[]
allReads = []
lenAllsites = []
 
i = 0

os.chdir(dir)
samfile = pysam.AlignmentFile(file, "rb")

for read in samfile.fetch(chrom, bedS,bedE):
    		
     		methylPos=methylpos(read.get_tag("XM"))
    		unmethylPos = unmethylpos(read.get_tag("XM"))
    		
    		sites = addLocus(methylPos, read.reference_start)
    		unsites = addLocus(unmethylPos, read.reference_start)
    		#print sites, "c not methylated:", unsites
    		newallSites = unsites + sites
    		n = len(newallSites)
	    	if n > max: max = n
	    	
	    	allReads = allReads + newallSites
	    	lenAllsites.append(len(newallSites))
	    	allReadsSites = allReadsSites + newallSites
	    	allReadsSites = set(allReadsSites)
	    	allReadsSites = list(allReadsSites)
samfile.close()
 

common_sitesNum = most_common(lenAllsites)

#print common_sitesNum
s = most_common(lenAllsites)
#print allReadsSites

top_nn = top_n(allReads, s)
top_nn.sort()

allReadsSites.sort()

j= 0
k = 0
name = re.sub("_R1_001_bismark_bt2.sorted.bam","", bam)
file2write = outputdir + '/'+ name + ".methylcounts.txt"
file3write = outputdir + '/'+ name + "_metrics.txt"
w = open(file2write,'w')
q = open(file3write, 'w')
string = 'H' + '\t'+ 'X'+ '\t'+ 'h' + '\t'+ 'x' + '\t'+ '\n'
str1 = 'SeqRead'+'\t'+ '\t'.join(chrom + str(e)  for e in top_nn)  + '\t'+ 'total' + '\n' 
w.write(str1)
q.write(string)
print chrom, bedS, bedE

samfile = pysam.AlignmentFile(bam, "rb" )
for read in samfile.fetch(chrom, bedS,bedE):

	vec =[]
	for jj in range(len(top_nn)):vec.append(0)
	methylPos=methylpos(read.get_tag("XM"))
	unmethylPos = unmethylpos(read.get_tag("XM"))
	allMethylC = methylPos + unmethylPos
	k = k + 1
	c = [ read.query_alignment_qualities[i] for i in allMethylC ]
	if  ( not(all(i >= 20 for i in c)) or read.mapping_quality < 20): continue
	j = j + 1
	methylsites = addLocus(methylPos, read.reference_start)
	unmethylsites = addLocus(unmethylPos, read.reference_start)
	HH = methylCitocineCHH(read.get_tag("XM"))
	XX = methylCitocineCHG(read.get_tag("XM"))
	x = unmethylCitocineCHH(read.get_tag("XM"))
	h = unmethylCitocineCHG(read.get_tag("XM"))
	a =str(len(HH))
	b = str(len(XX))
	c = str(len(h))
	d = str(len(x))
	
	meth = [top_nn.index(x) for x in methylsites if x in top_nn]
	unmeth = [top_nn.index(x) for x in unmethylsites if x in top_nn]
	for ii in meth:vec[ii]=1
	sumvec = sum(vec)
	vec2write = '\t'.join(str(e) for e in vec)
	total = str(sumvec)
	interval =  chrom + ':' + str(bedS)+'-'+str(bedE) 
	str3 =interval + '\t' +  vec2write+ '\t' + total + '\n'
	
	str4 =  a + '\t' + b +  '\t' +  c +  '\t' + d  + '\n'
	print str4
	print str3
	q.write(str4)
	w.write(str3)
samfile.close()

str5 = 'total reads = ' + str(k) + '\t' + 'high quality reads = ' + str(j) 
q.write(str5)
q.close()
w.close()



print j
print k
print top_nn










