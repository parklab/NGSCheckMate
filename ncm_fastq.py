import os
import sys
import math
import subprocess, time
import argparse
from argparse import RawTextHelpFormatter

global fastq1
global fastq2
global sub_rate
global desire_depth
global reference_length
global pattern_length
global maxthread
global PE
global bed_file
global outdir
global outfilename
global temp_out
global testsamplename

glob_scores = dict()    #Whole score
feature_list = dict()   #Each Feature List
label = []              #Samples
features = []           #dbSNP features
mean_depth = dict()
real_depth = dict()
sum_file = dict()
Family_flag = False

#Calculation of AVerages
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

#Calulation of Pearson Correlation
def pearson_def(x, y):
    assert len(x) == len(y)
    
    n = len(x)
    if n<20 :
        return 0
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
        
    # Remove devided by 0 cases     
    if math.sqrt(xdiff2 * ydiff2) ==0:
        return diffprod / (math.sqrt(xdiff2 * ydiff2) + 0.00001) 
        
        
    return diffprod / math.sqrt(xdiff2 * ydiff2)

# createDataSet
# base_dir : directory of files, bedFile: name of the bedFile
def createDataSetFromDir(base_dir, bedFile):
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if not file.endswith("ncm"):
                continue
                
#            if file.endswith("class_results.txt"):
#                continue
                
            link = root + '/' +  file
            f = open(link, "r")
  #          dbsnpf= open(bedFile,"r")
            depth = 0  
            count = 0
            real_count = 0 
 #           sum = 0
            
     #       file = file +"_" + order
            
            scores = dict()     # Scores of B-allel Frequencies
            #DBSNP ID collecting system
            for i in range(0,21039):
             #  temp = i.split('\t')
             #   ID = temp[0]
             #   scores[ID] = 0
                scores[str(i)] = 0
                count=count + 1
                    
            feature_list[file] = []
            #VCF file PROCESSING  and Generation of features  
            for line in f.readlines():
                if line.startswith("index"):
                    continue
                    
                temp = line.strip().split("\t")
                if temp[3] != "NA" and temp[3] != "vaf" and len(temp) > 3:
                    scores[temp[0]] = float(temp[3])
                    real = int(temp[1]) + int(temp[2])
                    depth = depth+ real
                    count = count + 1
                    if real > 0 :
                        real_count = real_count + 1   
                
                    feature_list[file].append(temp[0])
                    
            mean_depth[file] = depth / float(count) 
 #           print count
            real_depth[file] = depth / float(real_count)
 #           sum_file[file] = sum                        
                     
            for key in features:
                if glob_scores.has_key(file):
                    glob_scores[file].append(scores[key])
                else: 
                    glob_scores[file] = [scores[key]]
                          
 #           dbsnpf.close()
            f.close()            

    for key in sorted(glob_scores):
        label.append(key)    

# createDataSet
# base_dir : directory of files, bedFile: name of the bedFile
def createDataSetFromDir_test(base_dir, bedFile,order):
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if not file.endswith("ncm"):
                continue
                
#            if file.endswith("class_results.txt"):
#                continue
                
            link = root + '/' +  file
            f = open(link, "r")
  #          dbsnpf= open(bedFile,"r")
            depth = 0  
            count = 0
            real_count = 0 
 #           sum = 0
            
            file = file +"_" + order
            
            scores = dict()     # Scores of B-allel Frequencies
            #DBSNP ID collecting system
            for i in range(0,21039):
             #  temp = i.split('\t')
             #   ID = temp[0]
             #   scores[ID] = 0
                scores[str(i)] = 0
                count=count + 1
                    
            feature_list[file] = []
            #VCF file PROCESSING  and Generation of features  
            for line in f.readlines():
                if line.startswith("index"):
                    continue
                    
                temp = line.strip().split("\t")
                if temp[3] != "NA" and temp[3] != "vaf" and len(temp) > 3:
                    scores[temp[0]] = float(temp[3])
                    real = int(temp[1]) + int(temp[2])
                    depth = depth+ real
                    count = count + 1
                    if real > 0 :
                        real_count = real_count + 1   
                
                    feature_list[file].append(temp[0])
                    
            mean_depth[file] = depth / float(count) 
 #           print count
            real_depth[file] = depth / float(real_count)
 #           sum_file[file] = sum                        
                     
            for key in features:
                if glob_scores.has_key(file):
                    glob_scores[file].append(scores[key])
                else: 
                    glob_scores[file] = [scores[key]]
                          
 #           dbsnpf.close()
            f.close()            

    for key in sorted(glob_scores):
        label.append(key)    



def classifyNV(vec2Classify, p0Vec, p0S, p1Vec, p1S):    
    if abs(p0Vec - vec2Classify) + p0S > abs(p1Vec - vec2Classify) - p1S:
        return abs((abs(p0Vec - vec2Classify) +  p0S )/ (abs(p1Vec - vec2Classify) -  p1S )), 1
    else: 
        return abs((abs(p0Vec - vec2Classify) + p0S) / (abs(p1Vec - vec2Classify)  -  p1S)), 0        


def getPredefinedModel(depth):
     if Family_flag:
         if depth > 10:
             return 0.874546, 0.022211, 0.646256175, 0.021336239
         elif depth > 5:
             return 0.785249,0.021017, 0.598277053, 0.02253561
         elif depth > 2:
             return 0.650573, 0.018699,0.536020197, 0.020461932
         elif depth > 1:
             return 0.578386,0.018526, 0.49497342, 0.022346597
         elif depth > 0.5:
             return 0.529327,0.025785, 0.465275173, 0.028221203
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.529327,0.025785, 0.465275173, 0.028221203
     else:
         if depth > 10:
             return 0.874546, 0.022211, 0.310549, 0.060058
         elif depth > 5:
             return 0.785249,0.021017, 0.279778, 0.054104
         elif depth > 2:
             return 0.650573, 0.018699,0.238972, 0.047196
         elif depth > 1:
             return 0.578386,0.018526, 0.222322, 0.041186
         elif depth > 0.5:
             return 0.529327,0.025785, 0.217839, 0.040334
         else:
    #         print "Warning: Sample region depth is too low < 1"
             return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 30:
#         return 0.874546, 0.022211, 0.310549, 0.060058
#     elif depth > 10:
#         return 0.785249,0.021017, 0.279778, 0.054104
#     elif depth > 5:
#         return 0.650573, 0.018699,0.238972, 0.047196
#     elif depth > 2:
#         return 0.578386,0.018526, 0.222322, 0.041186
#     elif depth > 1:
#         return 0.529327,0.025785, 0.217839, 0.040334
#     else:
#         print "Warning: Sample region depth is too low < 1"
#         return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.1:
#        return 0.0351* depth + 0.5538, 0.02, 0.009977*depth + 0.216978, 0.045
#     else:
#        print "too low depth"
#        return 0.529327,0.025785, 0.217839, 0.040334
#     if depth > 0.5:
#        return 0.06315* (math.log(depth)) + 0.64903, 0.046154, 0.0005007*depth + 0.3311504,0.12216
#     else:
#        return 0.62036, 0.046154, 0.31785, 0.12216

def calAUC(predStrengths, classLabels):
    ySum = 0.0 #variable to calculate AUC
    cur = (1.0,1.0) #cursor
    numPosClas = sum(array(classLabels)==1.0)
    yStep = 1/float(numPosClas); xStep = 1/float(len(classLabels)-numPosClas)
    sortedIndicies = predStrengths.argsort()#get sorted index, it's reverse
    #loop through all the values, drawing a line segment at each point
    for index in sortedIndicies.tolist()[0]:
        if classLabels[index] == 1:
            delX = 0; delY = yStep;
        else:
            delX = xStep; delY = 0;
            ySum += cur[1]
        cur = (cur[0]-delX,cur[1]-delY)
    return ySum*xStep        

def plotROC(predStrengths, classLabels):
    import matplotlib.pyplot as plt
    cur = (1.0,1.0) #cursor
    ySum = 0.0 #variable to calculate AUC
    numPosClas = sum(array(classLabels)==1.0)
    yStep = 1/float(numPosClas); xStep = 1/float(len(classLabels)-numPosClas)
    sortedIndicies = predStrengths.argsort()#get sorted index, it's reverse
    fig = plt.figure()
    fig.clf()
    ax = plt.subplot(111)
    #loop through all the values, drawing a line segment at each point
    for index in sortedIndicies.tolist()[0]:
        if classLabels[index] == 1:
            delX = 0; delY = yStep;
        else:
            delX = xStep; delY = 0;
            ySum += cur[1]
        #draw line from cur to (cur[0]-delX,cur[1]-delY)
        ax.plot([cur[0],cur[0]-delX],[cur[1],cur[1]-delY], c='b')
        cur = (cur[0]-delX,cur[1]-delY)
    ax.plot([0,1],[0,1],'b--')
    plt.xlabel('False positive rate'); plt.ylabel('True positive rate')
    plt.title('ROC curves')
    ax.axis([0,1,0,1])
    plt.show()
    print "the Area Under the Curve is: ",ySum*xStep


def run_fastq_version():
    command = "./ngscheckmate_fastq "
    if sub_rate!= "":
            command = command + "-s " + sub_rate  + " "
    if desired_depth !="":
            command = command + "-d " + desired_depth + " "
    if reference_length !="":
            command = command + "-R " + reference_length + " "
    if pattern_length !="":
            command = command + "-L " + pattern_length + " "
    if maxthread !="":
            command = command + "-p " + maxthread + " "

    if PE == 1:
            command =  command  + "-1 "  + fastq1 + " -2 " + fastq2 +" " + bed_file +" > " + outdir + "/" + temp_out + ".ncm"

    print command

    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()



def classifying():
    AUCs =[]

    wholeFeatures = 50

    temp = []

    altFreqList = []
    keyList = []

    for key in sorted(glob_scores):
        altFreqList.append(glob_scores[key])
        keyList.append(key)

    dataSetSize = len(altFreqList)

    for i in range(0, dataSetSize):
        for j in range(0, dataSetSize):
            if i!=j:
                temp.append([keyList[i],keyList[j]])

    for iterations in range(49,wholeFeatures):

        samples = []
        numFeatures = iterations

        count = 0

        for i in range(0,len(temp)):
            tempA = set(feature_list[temp[i][0].strip()])
            tempB = set(feature_list[temp[i][1].strip()])

            selected_feature = tempA.intersection(tempB)
            
            vecA = []
            vecB = []
            
            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(glob_scores[temp[i][0].strip()][idx])
                    vecB.append(glob_scores[temp[i][1].strip()][idx])
                idx = idx + 1
            
            
            
            distance = pearson_def(vecA, vecB)
            samples.append(distance)
            
    predStrength = []
    training_flag =0
####0715 Append

    output_matrix_f = open(outdir + "/output_corr_matrix.txt","w")
    output_matrix = dict()
    
    if out_tag!="stdout":
    	out_f = open(outdir + "/" + out_tag + ".txt","w")

    for i in range(0, len(samples)):
        output_matrix[temp[i][0]] = dict()
        for j in range(0,len(samples)):
            output_matrix[temp[i][0]][temp[j][0]] = 0

    if training_flag == 1:
        #make training set
        for i in range(0,len(samples)):
            trainMatrix= []
            trainCategory = []
            for j in range(0, len(samples)):
                if i==j:
                    continue
                else:
                    trainMatrix.append(samples[j])
                    trainCategory.append(classLabel[j])
            #training samples in temp
            #p0V, p1V, pAb = trainNB0(array(trainMatrix),array(trainCategory))
            p1V,p1S, p0V, p0S = trainNV(array(trainMatrix),array(trainCategory))
            result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
            if result[1] == 1:
                print str(temp[i][0]) + '\tsample is matched to\t',str(temp[i][1]),'\t', samples[i]
            predStrength.append(result[0])
#            AUCs.append(calAUC(mat(predStrength),classLabel))
#            plotROC(mat(predStrength),classLabel)
#            print AUCs
    else :
        for i in range(0,len(samples)):
            depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])
            p1V,p1S, p0V, p0S = getPredefinedModel(depth)
            result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
            if result[1] ==1:
                output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                if out_tag=="stdout":
                    print str(temp[i][0][:-4]) + '\tmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                else :
                    out_f.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
            #print sum_file[temp[i][0]],sum_file[temp[i][1].strip()]
            predStrength.append(result[0])
#            AUCs.append(calAUC(mat(predStrength),classLabel))
#            plotROC(mat(predStrength),classLabel)
#            print AUCs
        #testing sample is samples
    output_matrix_f.write("sample_ID")
    for key in output_matrix.keys():
        output_matrix_f.write("\t" + key[0:key.index('.')])
    output_matrix_f.write("\n")

    for key in output_matrix.keys():
        output_matrix_f.write(key[0:key.index('.')])
        for otherkey in output_matrix.keys():
            output_matrix_f.write("\t" + str(output_matrix[key][otherkey]))
        output_matrix_f.write("\n")   
        
    output_matrix_f.close()         
    if out_tag!="stdout":
    	out_f.close()   


def classifying_test():
    AUCs =[]

    wholeFeatures = 50

    temp = []

    keyF = open(samplefilename,'r')
    temp =[]

    for k in outF.readlines():
        keyfile = k.split(":")
        keyfile[0] = keyfile[0].strip() + "_1"
        keyfile[1] = keyfile[1].strip() + "_2"
        temp.append(keyfile)
    keyF.close()

    for iterations in range(49,wholeFeatures):

        samples = []
        numFeatures = iterations

        count = 0

        for i in range(0,len(temp)):
            tempA = set(feature_list[temp[i][0].strip()])
            tempB = set(feature_list[temp[i][1].strip()])

            selected_feature = tempA.intersection(tempB)
            
            vecA = []
            vecB = []
            
            idx = 0
            for k in features:
                if k in selected_feature:
                    vecA.append(glob_scores[temp[i][0].strip()][idx])
                    vecB.append(glob_scores[temp[i][1].strip()][idx])
                idx = idx + 1
            
            distance = pearson_def(vecA, vecB)
            samples.append(distance)
            
    predStrength = []
    training_flag =0
####0715 Append

    output_matrix_f = open(outdir + "/output_corr_matrix.txt","w")
    output_matrix = dict()
    
    if out_tag!="stdout":
        out_f = open(outdir + "/" + out_tag + ".txt","w")

    for i in range(0, len(samples)):
        output_matrix[temp[i][0]] = dict()
        for j in range(0,len(samples)):
            output_matrix[temp[i][0]][temp[j][0]] = 0

    if training_flag == 1:
        #make training set
        for i in range(0,len(samples)):
            trainMatrix= []
            trainCategory = []
            for j in range(0, len(samples)):
                if i==j:
                    continue
                else:
                    trainMatrix.append(samples[j])
                    trainCategory.append(classLabel[j])
            #training samples in temp
            #p0V, p1V, pAb = trainNB0(array(trainMatrix),array(trainCategory))
            p1V,p1S, p0V, p0S = trainNV(array(trainMatrix),array(trainCategory))
            result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
            if result[1] == 1:
                print str(temp[i][0]) + '\tsample is matched to\t',str(temp[i][1]),'\t', samples[i]
            predStrength.append(result[0])
#            AUCs.append(calAUC(mat(predStrength),classLabel))
#            plotROC(mat(predStrength),classLabel)
#            print AUCs
    else :
        for i in range(0,len(samples)):
            depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])
            p1V,p1S, p0V, p0S = getPredefinedModel(depth)
            result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
            if result[1] ==1:
                output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                if out_tag=="stdout":
                    print str(temp[i][0][:-4]) + '\tmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                else :
                    out_f.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
            #print sum_file[temp[i][0]],sum_file[temp[i][1].strip()]
            predStrength.append(result[0])
#            AUCs.append(calAUC(mat(predStrength),classLabel))
#            plotROC(mat(predStrength),classLabel)
#            print AUCs
        #testing sample is samples
    output_matrix_f.write("sample_ID")
    for key in output_matrix.keys():
        output_matrix_f.write("\t" + key[0:key.index('.')])
    output_matrix_f.write("\n")

    for key in output_matrix.keys():
        output_matrix_f.write(key[0:key.index('.')])
        for otherkey in output_matrix.keys():
            output_matrix_f.write("\t" + str(output_matrix[key][otherkey]))
        output_matrix_f.write("\n")   
        
    output_matrix_f.close()         
    if out_tag!="stdout":
        out_f.close()   

def generate_R_scripts():
    r_file = open(outdir + "/r_script.r","w")
    cmd = "output_corr_matrix <- read.delim(\"" + outdir +  "/output_corr_matrix.txt\")\n"
    cmd = cmd + "data = output_corr_matrix\n"
    cmd = cmd + "d3 <- as.dist((1 - data[,-1]))\n"
    cmd = cmd + "clust3 <- hclust(d3, method = \"average\")\n"
    cmd = cmd + "pdf(\"" +outdir+ "/" + pdf_tag + ".pdf\", width=10, height=7)\n"
    cmd = cmd + "op = par(bg = \"gray85\")\n"
    cmd = cmd + "par(plt=c(0.05, 0.95, 0.5, 0.9))\n"
    cmd = cmd + "plot(clust3, lwd = 2, lty = 1,cex=0.8, xlab=\"Samples\", sub = \"\",  ylab=\"Distance (1-Pearson correlation)\",hang = -1, axes = FALSE)\n"
    cmd = cmd + "axis(side = 2, at = seq(0, 1, 0.2), labels = FALSE, lwd = 2)\n"
    cmd = cmd + "mtext(seq(0, 1, 0.2), side = 2, at = seq(0, 1, 0.2), line = 1,   las = 2)\n"
    cmd = cmd + "dev.off()\n"
    r_file.write(cmd)
    r_file.close()


def run_R_scripts():
    command = "R CMD BATCH " + outdir + "/r_script.r"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()

def remove_internal_files():
    if outdir.find("*"):
        sys.exit()


    command = "rm -rf " + outdir + "/output_corr_matrix.txt"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()
    command = "rm -rf " + outdir + "/r_script.r"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()
    command = "rm -rf " + outdir + "/r_script.r.Rout"
    proc = subprocess.Popen(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code = proc.wait()

if __name__ == '__main__':
    sub_rate = ""
    desired_depth = ""
    reference_length =""
    pattern_length = ""
    maxthread =""
    PE = 0
    fastq1 = ""
    fastq2 = ""
    testsamplename = ""

    help = """
    Ensuring Sample Identity v0.8
    Usage : python ./ngscheckmate_fastq <options> -1 fastqfile1 [-2 fastqfile2]  patternfile

        Input arguments (required)
          patternfile : a text file with sequences flanking representative snv sites, along with markers indicating the snv index and whether the sequence represents reference or alternative allele.
          fastqfile1 : see below 'Options'.

        Options
          -s, --ss <subsampling_rate> : subsampling rate (default 1.0)
          -d, --depth <desired_depth> : as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data.
          -R, --reference_length <reference_length> : The reference length (default : 3E9) to be used for computing subsampling rate. If the data is NOT WGS from human, and if you're using the -d option, it is highly recommended to specify the reference length. For instance, if your data is human RNA-seq, the total reference length could be about 3% of the human genome, which can be set as 1E8.
          -L, --pattern_length <pattern_length> : The length of the flanking sequences being used to identify SNV sites. Default is 21bp. It is recommended not to change this value, unless you have created your own pattern file with a different pattern length.
          -p, --maxthread <number_of_threads> : number of threads to use (default : 1 )

          -O, --output_dir <out_dir> : A directory of output files. <default : ./>
          -N, --Name <result_name> : a prefix of name of result file. <default : output>
          -I, --input_file <absolute_path_fastq_files_list> : List of Absolute paths of fastq files (single file per line, if second fastq file is exist, write fastq files using tab delimited. (required)
                                                             ex) fastq_file_R1  fastq_file_R2   sample_ID
          -f, --family-considered mode : Reference correlations of family considered models
          -c, --answerset <answer_file> : Input answer file which have pair list of file names user want to examine (default : all possible combinations)


    Sejoon Lee, Soo Lee, Eunjung Lee, 2015
            """

    parser = argparse.ArgumentParser(description=help, formatter_class=RawTextHelpFormatter)

#    group_type = parser.add_mutually_exclusive_group(required=True)
#    group_type.add_argument()
#    group = parser.add_mutually_exclusive_group(required=True)
#    group.add_argument('-v','--vcf',metavar='VCF_list',dest='vcf_files_list',action='store', help='VCF files from samtools mpileup and bcftools')
#    group.add_argument('-d','--dir',metavar='VCF_dir',dest='vcf_files_dir',action='store', help='VCF files from samtools mpileup and bcftools')

    parser.add_argument('-f','--family_cutoff',dest='family_cutoff',action='store_true', help='apply strict correlation threshold to remove family cases')
    parser.add_argument('-bed','--bed',metavar='feature bed file',required=True,dest='bed_file',action='store', help='bed file')
    parser.add_argument('-s','--ss',metavar='subsampling_rate',dest='sub_rate',action='store', help='subsampling rate (default 1.0)')
    parser.add_argument('-d','--depth',metavar='desired_depth',dest='desired_depth',action='store', help='as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data')
    parser.add_argument('-R','--reference_length',metavar='reference_length',dest='reference_length',action='store', help="The reference length (default : 3E9) to be used for computing subsampling rate. If the data is NOT WGS from human, and if you aree using the -d option, it is highly recommended to specify the reference length. For instance, if your data is human RNA-seq, the total reference length could be about 3percent of the human genome, which can be set as 1E8.")
    parser.add_argument('-L','--pattern_length',metavar='pattern_length',dest='pattern_length',action='store', help='The length of the flanking sequences being used to identify SNV sites. Default is 21bp. It is recommended not to change this value, unless you have created your own pattern file with a different pattern length.')
    parser.add_argument('-p','--maxthread',metavar='maxthread',dest='maxthread',action='store', help='number of threads to use (default : 1 )')

    parser.add_argument('-O','--outdir',metavar='output_dir',dest='outdir',action='store', help='directory name for temp and output files')
    parser.add_argument('-N','--outfilename',metavar='output_filename',dest='outfilename',action='store',default="output",help='OutputFileName ( default : output ), -N filename')
    parser.add_argument('-I','--input',metavar='input_file_name',required=True,dest='inputfilename',action='store',help='Inputfile name that contains fastq file names, -I filename')

    parser.add_argument('-t','--testsamplename',metavar='test_samplename',dest='testsamplename',action='store',help='file including test sample namses  with ":" delimeter (default : all combinations of samples), -t filename')


    args=parser.parse_args()

    bed_file = args.bed_file
    outdir = args.outdir
    outfilename = args.outfilename

    if args.sub_rate != None:
        sub_rate = args.sub_rate
    if args.desired_depth != None:
        desired_depth = args.desired_depth
    if args.reference_length != None:
        reference_length = args.reference_length
    if args.pattern_length != None:
        pattern_length = args.pattern_length
    if args.maxthread != None:
        maxthread = args.maxthread

    if args.family_cutoff:
        Family_flag=True

    with open(args.inputfilename,'r') as F:
        for line in F.xreadlines():
                temp = line.strip().split("\t")

                if len(temp) == 3:
                        PE = 1
                        fastq1 = temp[0]
                        fastq2 = temp[1]
                        temp_out = temp[2]
                        run_fastq_version()
                elif len(temp) == 2:
                        PE = 0
                        fastq1 = temp[0]
                        temp_out = temp[1]
                        run_fastq_version()
                else:
                        print "Input File Error: Each line should be contain one or two fastq files name with tab delimited"
                        print line.strip()
                        print "upper format is invalid"

    # set directories
    base_dir = outdir

    #base_dir = "/data/users/sjlee/valid_qc/WGS/SNP/MATCH/"

    #bedFile = "/data/users/sjlee/qc/disctinct_9755.bed"
    bedFile = bed_file
    #outFileName = "/data/users/sjlee/valid_qc/WGS/SNP/MATCH_CLASS/wgs_CL.txt"
#    outFileName = args.class_file
    out_tag = outfilename
#    key_feature_F = "/data/users/sjlee/qc/vcf_generator/feature_selection/Distinct_9755_features.txt"

    
#    outCL = open(outFileName[:outFileName.index('.')]+'.class','r')

#    classLabel=[]
#    for i in outCL.readlines():
#        classLabel.append(int(i.strip()))

    #key_order = open(key_feature_F,'r')
    key_order = open(bedFile,'r')

    fastq = 1
    
    if fastq == 0:
        for i in key_order.readlines():
            temp = i.split('\t')
            features.append(str(temp[0])+"_"+str(temp[2]))
            
    if fastq == 1:
        for i in range(0,21039):
            features.append(str(i))     

    if args.testsamplename != None:
        testsamplename = args.testsamplename
        print "Generate Data Set from " + outdir + "\nusing this bed file : " + bedFile
        createDataSetFromDir_test(outdir,bedFile,"1")
        createDataSetFromDir_test(outdir,bedFile,"2")
        classifying_test()
    else:
        print "Generate Data Set from " + outdir + "\nusing this bed file : " + bedFile
        createDataSetFromDir(outdir,bedFile)
        classifying()


#    print "Generate Data Set from " + outdir + "\nusing this bed file : " + bedFile
#    createDataSetFromList(outdir,bedFile)

#    if args.method == "clustering":
#        print "Classifying data set based on kNN ",str(args.KNN)
#        clustering(int(args.KNN))
#    elif args.method =="classifying":
    

  
#  if args.PDF_flag != None:
    pdf_tag = outfilename
    generate_R_scripts()
    run_R_scripts()
#   remove_internal_files()