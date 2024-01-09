#!/usr/bin/env python
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
global nodeptherror
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
Nonzero_flag = False

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
            if float(real_count) == 0:
                real_depth[file] = depth / float(count)
            else:
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
            if float(real_count) == 0:
                real_depth[file] = depth / float(count)
            else:
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
    if abs(p0Vec - vec2Classify) - p0S > abs(p1Vec - vec2Classify) - p1S:
        return abs((abs(p0Vec - vec2Classify) - p0S )/ (abs(p1Vec - vec2Classify) -  p1S )), 1
    else: 
        return abs((abs(p0Vec - vec2Classify) - p0S) / (abs(p1Vec - vec2Classify)  -  p1S)), 0        


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
    INSTALL_DIR=""
    if "NCM_HOME" in os.environ.keys():
        INSTALL_DIR=os.environ['NCM_HOME'] + "/"
    else :
        print "WARNNING : NCM_HOME is not defined yet. Therefore, program will try to search ngscheckmate_fastq file from the current directory"
        INSTALL_DIR="./"

    command = INSTALL_DIR + "ngscheckmate_fastq "
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
    if nodeptherror !="":
            command = command + "-j " + nodeptherror + " "

    if PE == 1:
            command =  command  + "-1 "  + fastq1 + " -2 " + fastq2 +" " + bed_file +" > " + outdir + "/" + temp_out + ".ncm"
    if PE == 0:
	        command = command + "-1 " + fastq1  +" " + bed_file +" > " + outdir + "/" + temp_out + ".ncm"

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

    filter_list = []

    for i in range(0, dataSetSize):
        for j in range(0, dataSetSize):
            if i!=j:
                if keyList[j] not in filter_list:
                    temp.append([keyList[i],keyList[j]])
        filter_list.append(keyList[i])

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
            out_f = open(outdir + "/" + out_tag + "_all.txt","w")
            out_matched = open(outdir + "/" + out_tag + "_matched.txt","w")

        for i in range(0, len(keyList)):
            output_matrix[keyList[i]] = dict()
            for j in range(0,len(keyList)):
                output_matrix[keyList[i]][keyList[j]] = 0

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
                depth = 0 
                if Nonzero_flag:
                    depth = min(real_depth[temp[i][0].strip()],real_depth[temp[i][1].strip()])
                else:
                    depth = min(mean_depth[temp[i][0].strip()],mean_depth[temp[i][1].strip()])
  
                p1V,p1S, p0V, p0S = getPredefinedModel(depth)
                result = classifyNV(samples[i],p0V,p0S, p1V, p1S)
                if result[1] ==1:
                    output_matrix[temp[i][0].strip()][temp[i][1].strip()] = samples[i]
                    output_matrix[temp[i][1].strip()][temp[i][0].strip()] = samples[i]
                    if out_tag=="stdout":
                        print str(temp[i][0][:-4]) + '\tmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                        out_matched.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')               
                else:
                    if out_tag=="stdout":
                        print str(temp[i][0][:-4]) + '\tunmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                    else :
                        out_f.write(str(temp[i][0][:-4]) + '\tunmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
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
            out_matched.close()   


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
        out_f = open(outdir + "/" + out_tag + "_all.txt","w")
        out_matched = open(outdir + "/" + out_tag + "_matched.txt","w")

    for i in range(0, len(keyList)):
        output_matrix[keyList[i]] = dict()
        for j in range(0,len(keyList)):
            output_matrix[keyList[i]][keyList[j]] = 0

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
                output_matrix[temp[i][1].strip()][temp[i][0].strip()] = samples[i]
                if out_tag=="stdout":
                    print str(temp[i][0][:-4]) + '\tmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                else :
                    out_f.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
                    out_matched.write(str(temp[i][0][:-4]) + '\tmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')               
            else:
                if out_tag=="stdout":
                    print str(temp[i][0][:-4]) + '\tunmatched\t',str(temp[i][1][:-4]),'\t', round(samples[i],4),'\t',round(depth,2)
                else :
                    out_f.write(str(temp[i][0][:-4]) + '\tunmatched\t' + str(temp[i][1][:-4])  + '\t'+  str(round(samples[i],4)) + '\t' + str(round(depth,2)) + '\n')
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
        out_matched.close()   

def generate_R_scripts():
    r_file = open(outdir + "/r_script.r","w")
    if len(feature_list)==0:
       r_file.close()
    else :
       cmd = "output_corr_matrix <- read.delim(\"" + outdir +  "/output_corr_matrix.txt\")\n"
       cmd = cmd + "data = output_corr_matrix\n"
       cmd = cmd + "d3 <- as.dist((1 - data[,-1]))\n"
       cmd = cmd + "clust3 <- hclust(d3, method = \"average\")\n"
       if len(feature_list) < 5:
           cmd = cmd + "pdf(\"" +outdir+ "/" + pdf_tag + ".pdf\", width=10, height=7)\n"
       else:
           cmd = cmd + "pdf(\"" +outdir+ "/" + pdf_tag + ".pdf\", width="+str(math.log10(len(feature_list))*10) +", height=7)\n"
       cmd = cmd + "op = par(bg = \"gray85\")\n"
       cmd = cmd + "par(plt=c(0.05, 0.95, 0.2, 0.9))\n"
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

def output_filter():
    success_set_M = []
    success_set_U = []
    failure_set_M = []
    failure_set_U = []

    with open(outdir + "/" + out_tag + "_all.txt","r") as F:
        for line in F.readlines():
            temp = line.strip().split('\t')
            
            sample1 = temp[0]
            sample2 = temp[2]
            
            match = temp[1]
            
            if match == "matched":
                if sample1[sample1.index("TCGA"):sample1.index("TCGA")+12] == sample2[sample2.index("TCGA"):sample2.index("TCGA")+12] :
                    success_set_M.append(line)
                else:
                    failure_set_M.append(line)
            elif match == "unmatched":
                if sample1[sample1.index("TCGA"):sample1.index("TCGA")+12] == sample2[sample2.index("TCGA"):sample2.index("TCGA")+12] :
                    failure_set_U.append(line)
                else:
                    success_set_U.append(line)        
              
    Matched_file = open(outdir + "/" + out_tag + "_matched.txt",'w') 

    for i in success_set_M:
        Matched_file.write(i)
    for i in failure_set_M:
        Matched_file.write(i)  
    
    Matched_file.close()

    problem_file = open(outdir + "/" + out_tag + "_problematic.txt",'w')

    for i in failure_set_M:
        problem_file.write(i)
    for i in failure_set_U:
        problem_file.write(i)

    problem_file.close()

    Summary_file = open(outdir + "/" + out_tag + "_summary.txt",'w')
    
 

    ## paired cluster - only failed things
    Summary_file.write("###########################################\n")
    Summary_file.write("###  Problematic clusters of same orgins ##\n")
    Summary_file.write("###########################################\n\n")

    cluster = dict()

    result_set = failure_set_M + success_set_M

    for line in result_set:
        temp = line.strip().split('\t')
        flag = 0
        for key in cluster:
            if temp[0] in cluster[key]:
                cluster[key].add(temp[2])
                flag = 1
                break
            elif temp[2] in cluster[key]:
                cluster[key].add(temp[0])
                flag = 1
                break
        
        if flag == 0:
            cluster[temp[0]] = set()
            cluster[temp[0]].add(temp[0])
            cluster[temp[0]].add(temp[2])
            
            
    count = 0 
    for key in cluster:
        temp_list = []
        flag = 0
        for data in cluster[key]:
            temp_list.append(data)
            sample1 = temp_list[0]
            ID = sample1[sample1.index("TCGA"):sample1.index("TCGA")+12]
            
            for sample1 in cluster[key]:
                if ID != sample1[sample1.index("TCGA"):sample1.index("TCGA")+12]:
                    flag = 1

              

        if flag == 1:
            count = count + 1
            Summary_file.write("Cluster " + str(count) + "\n")
              
            for data in cluster[key]:
                Summary_file.write(data + "\n")
            Summary_file.write("\n")

                
    ## Singleton
    Summary_file.write("\n")
    Summary_file.write("###########################################\n")
    Summary_file.write("############### Singleton #################\n")
    Summary_file.write("###########################################\n\n")

    final_set = set()
    filter_set = set()

    result_set = failure_set_U

    for line in result_set:
        temp = line.strip().split('\t')
        
        final_set.add(temp[0])
        final_set.add(temp[2])
        
        flag = 0
        for key in cluster:
            if temp[0] in cluster[key]:
                filter_set.add(temp[0])
            elif temp[2] in cluster[key]:
                filter_set.add(temp[2])
                


    for i in final_set.difference(filter_set):
        Summary_file.write(i + "\n")

    Summary_file.close()



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
    nodeptherror = ""

    help = """
    NGSCheckMate v1.0
    Usage : python ncm_fastq.py -l INPUT_LIST_FILE -pt PT_FILE -O OUTPUT_DIR [options]
            python ncm_fastq.py -l FASTQ_list.txt -pt ./SNP/SNP.pt -O ./ncm_fastq_output -p 4 -f
            python ncm_fastq.py -l FASTQ_list.txt -pt ./SNP/SNP.pt -O ./ncm_fastq_output -p 4 -f -nz

        Input arguments (required)
          -l  FILE      A text file that lists input fastq (or fastq.gz) files and sample names (one per line; see Input file format)
          -pt FILE      A binary pattern file (.pt) that lists flanking sequences of selected SNPs (included in the package; SNP/SNP.pt)
          -O  DIR       An output directory

        Options
          -N PREFIX     A prefix for output files (default: "output")

          -f            Use strict VAF correlation cutoffs. Recommended when your data may include
                        related individuals (parents-child, siblings)

          -nz           Use the mean of non-zero depths across the SNPs as a reference depth
                        (default: Use the mean depth across all the SNPs)

          -s FLOAT      The read subsampling rate (default: 1.0)

          -d INT        The target depth for read subsampling. NGSCheckMate calculates a subsampling rate based on this target depth.

          -R INT        The length of the genomic region with read mapping (default: 3E9) used to compute subsampling rate.
                        If your data is NOT human WGS and you use the -d option,
                        it is highly recommended that you specify this value.

          -L INT        The length of the flanking sequences of the SNPs (default: 21bp).
                        It is not recommended that you change this value unless you create your own pattern file (.pt) with a different length.
                        See Supporting Scripts for how to generate your own pattern file.

          -p INT        The number of threads (default: 1)
            """

    parser = argparse.ArgumentParser(description=help, formatter_class=RawTextHelpFormatter)

#    group_type = parser.add_mutually_exclusive_group(required=True)
#    group_type.add_argument()
#    group = parser.add_mutually_exclusive_group(required=True)
#    group.add_argument('-v','--vcf',metavar='VCF_list',dest='vcf_files_list',action='store', help='VCF files from samtools mpileup and bcftools')
#    group.add_argument('-d','--dir',metavar='VCF_dir',dest='vcf_files_dir',action='store', help='VCF files from samtools mpileup and bcftools')

    parser.add_argument('-f','--family_cutoff',dest='family_cutoff',action='store_true', help='apply strict correlation threshold to remove family cases')
    parser.add_argument('-pt','--pt',metavar='feature pattern file',required=True,dest='bed_file',action='store', help='pattern file')
    parser.add_argument('-s','--ss',metavar='subsampling_rate',dest='sub_rate',action='store', help='subsampling rate (default 1.0)')
    parser.add_argument('-d','--depth',metavar='desired_depth',dest='desired_depth',action='store', help='as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data')
    parser.add_argument('-R','--reference_length',metavar='reference_length',dest='reference_length',action='store', help="The reference length (default : 3E9) to be used for computing subsampling rate.")
    parser.add_argument('-L','--pattern_length',metavar='pattern_length',dest='pattern_length',action='store', help='The length of the flanking sequences being used to identify SNV sites. Default is 21bp.\nIt is recommended not to change this value, unless you have created your own pattern file with a different pattern length.')
    parser.add_argument('-p','--maxthread',metavar='maxthread',dest='maxthread',action='store', help='number of threads to use (default : 1 )')
    parser.add_argument('-j','--nodeptherror',metavar='nodeptherror',dest='nodeptherror',action='store', help='in case estimated subsampling rate is larger than 1, do not stop but reset it to 1 and continue')

    parser.add_argument('-O','--outdir',metavar='output_dir',dest='outdir',action='store', help='directory name for temp and output files')
    parser.add_argument('-N','--outfilename',metavar='output_filename',dest='outfilename',action='store',default="output",help='OutputFileName ( default : output ), -N filename')
    parser.add_argument('-l','--list',metavar='input_file_list',required=True,dest='inputfilename',action='store',help='Inputfile name that contains fastq file names, -I filename')
    parser.add_argument('-nz','--nonzero',dest='nonzero_read',action='store_true',help='Use non-zero mean depth of target loci as reference correlation. (default: Use mean depth of all target loci)')

    parser.add_argument('-t','--testsamplename',metavar='test_samplename',dest='testsamplename',action='store',help='file including test sample namses  with ":" delimeter (default : all combinations of samples), -t filename.\n-t option is for the previous NGSCheckMate version. No longer used.')


    args=parser.parse_args()

    bed_file = args.bed_file
    outdir = args.outdir
    outfilename = args.outfilename

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

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
    if args.nodeptherror != None:
        nodeptherror = args.nodeptherror

    if args.family_cutoff:
        Family_flag=True
    if args.nonzero_read:
        Nonzero_flag=True

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
#    output_filter()
    pdf_tag = outfilename
    generate_R_scripts()
    run_R_scripts()
#   remove_internal_files()
