# -*- coding: utf-8 -*-
# @Author  : LiQun
# @Email   : liqun95@163.com
# @Time    : 2021/1/19 19:33
# @File    : extractDenovo.py
# @Software: PyCharm

######################
# 需要以下信息
# family information
# 'bcftools query -l file.vcf'
# FID:familyID	IID:IndividualID	RID:relativeID  LOC:location
# 家系三人/两人+亲属
######################

######################
import os
import argparse
parser = argparse.ArgumentParser(description='extract de novo')
parser.add_argument('--vcf', type=str, required=True)
parser.add_argument('--ped', type=str, required=True)
parser.add_argument('--out', type=str, required=True)
args = parser.parse_args()
######################
sampleInfoFile = args.ped
inVcfile = args.vcf
outDir = args.out
folder = os.path.exists(outDir)
if not folder:
    os.makedirs(outDir)
######################
# 读取数据
# sampleInfo = open('samples.txt','r').readlines()
# inVcf = open('raw.vcf','r').readlines()
sampleInfo = open(sampleInfoFile,'r').readlines()
inVcf = open(inVcfile,'r').readlines()
######################

######################
# 函数式
# 提取基因型信息GT
def extractGT(str):
    geneType = str.split(':')[0]
    return geneType

# 提取测序深度信息
def extractDP(str):
    readDepth = str.split(':')[2]
    if readDepth == ".":
        readDepth = 0
    return readDepth

# 根据家系ID提取位置信息，包括父母与患者
def extractFamily(fid):
    sampleInfoFile = sampleInfo
    loc = {"patientLoc":0,"maternalLoc":0,"paternalLoc":0}
    for line in sampleInfoFile:
        line_sets = line.split('\t')
        if line_sets[0].strip() == str(fid):
            if line_sets[2] == "1":
                loc["patientLoc"] = line_sets[3].strip()
            if line_sets[2] == "2":
                loc["maternalLoc"] = line_sets[3].strip()
            if line_sets[2] == "3":
                loc["paternalLoc"] = line_sets[3].strip()
    return loc['patientLoc'],loc['maternalLoc'],loc['paternalLoc']

# 根据家系ID提取Format信息
def extractFormat(fid,linevcf):
    linevcf_sets = linevcf.split('\t')
    patientLoc1,maternalLoc1,paternalLoc1 = extractFamily(fid)
    #print(patientLoc1,maternalLoc1,paternalLoc1)
    patientLoc = int(patientLoc1) + 8
    if int(maternalLoc1) == 0:
        maternalLoc = 0
        maternalFormat = ".:.:.:.:."
    else:
        maternalLoc = int(maternalLoc1) + 8
        maternalFormat = linevcf_sets[maternalLoc]
    if int(paternalLoc1) == 0:
        paternalLoc = 0
        paternalFormat = ".:.:.:.:."
    else:
        paternalLoc = int(paternalLoc1) + 8
        paternalFormat = linevcf_sets[paternalLoc]
    #print(patientLoc,maternalLoc,paternalLoc)
    patientFormat = linevcf_sets[patientLoc]
    #maternalFormat = linevcf_sets[maternalLoc]
    #paternalFormat = linevcf_sets[paternalLoc]
    #print(patientFormat+'\t'+maternalFormat+'\t'+paternalFormat)
    return extractGT(patientFormat),extractGT(maternalFormat),extractGT(paternalFormat),patientLoc,maternalLoc,paternalLoc,extractDP(patientFormat),extractDP(maternalFormat),extractDP(paternalFormat)

# 根据家系ID提取亲属信息
def extractRelative(fid):
    tmp = [1,2,3]
    relative = []
    sampleInfoFile = sampleInfo
    for line in sampleInfoFile:
        line_sets = line.split('\t')
        if line_sets[0].strip() == str(fid):
            if int(line_sets[2]) not in tmp:
                relative.append(int(line_sets[3]))
    return relative

# 根据家系ID提取样本ID
def extractFamilyID(fid):
    tmp = []
    familysample = {}
    familysampleID = []
    sampleInfoFile = sampleInfo
    for line in sampleInfoFile:
        line_sets = line.split('\t')
        if line_sets[0].strip() == str(fid):
            if int(line_sets[2]) not in tmp:
                familysample[line_sets[2]] = line_sets[1]
                tmp.append(int(line_sets[2]))
    sortKeys = sorted(familysample.keys())
    for i in sortKeys:
        familysampleID.append(familysample[str(i)])
    return familysampleID
######################

######################
# 提取De novo位点
# 存储家系信息
familyID=[]
for line in sampleInfo:
    line_sets = line.split('\t')
    fid = line_sets[0].strip()
    if fid not in familyID:
        familyID.append(fid)
familyID.remove('FID')
######################

# 添加输出文件的头文件
for fid in familyID:
    fileName = outDir + "/out_" + fid + ".vcf"
    originHeader = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    relativeID = extractFamilyID(fid)
    num=len(relativeID)
    relativeInfo=relativeID[0]
    for i in range(1,num):
        relativeInfo+="\t"+relativeID[i].strip()
    header = originHeader + "\t" + relativeInfo
    with open(fileName,'a') as f:
        f.writelines(header+'\n')

######################
# main
# 按照家系进行遍历
for fid in familyID:
    print("family:"+fid)
    for line in inVcf:
        line_sets = line.split("\t")
        if line[0] == "#":
            pass
        else:
            if line.find("PASS") != -1:
                patientGT,maternalGT,paternalGT,patientLOC,maternalLOC,paternalLOC,patientDP,maternalDP,paternalDP=extractFormat(fid,line)
                if int(maternalLOC) == 0:
                    #print(line)
                    if int(patientDP) >= 10 and int(paternalDP) >=10:
                        if (patientGT == "0/1" or patientGT == "0|1" or patientGT == "1/1" or patientGT == "1|1" or patientGT == "1/0" or patientGT == "1|0") and (paternalGT == "0/0" or paternalGT == "0|0"):
                            sql = line_sets[0].strip() + "\t" + line_sets[1].strip() + "\t" + line_sets[2].strip() + "\t" + line_sets[3].strip()
                            sql = sql + "\t" + line_sets[4].strip() + "\t" + line_sets[5].strip() + "\t" + line_sets[6].strip() + "\t" + line_sets[7].strip()
                            sql = sql + "\t" + line_sets[8].strip()
                            sql = sql + "\t" + line_sets[int(patientLOC)].strip() + "\t" + line_sets[int(paternalLOC)].strip()
                            fileName = outDir + "/out_" + fid + ".vcf"
                            relativeLOC = extractRelative(fid)
                            if len(relativeLOC) != 0:
                                sql1 = line_sets[int(relativeLOC[0] + 8)].strip()
                                for i in range(1, len(relativeLOC)):
                                    sql1 = sql1 + "\t" + line_sets[int(relativeLOC[i] + 8)].strip()
                                sql = sql + "\t" + sql1
                            with open(fileName, "a") as f:
                                f.writelines(sql + "\n")
                        else:
                            pass
                if int(paternalLOC) == 0:
                    #print(line)
                    if int(patientDP) >= 10 and int(maternalDP) >= 10:
                        if (patientGT == "0/1" or patientGT == "0|1" or patientGT == "1/1" or patientGT == "1|1" or patientGT == "1/0" or patientGT == "1|0") and (maternalGT == "0/0" or maternalGT == "0|0"):
                            sql = line_sets[0].strip() + "\t" + line_sets[1].strip() + "\t" + line_sets[2].strip() + "\t" + line_sets[3].strip()
                            sql = sql + "\t" + line_sets[4].strip() + "\t" + line_sets[5].strip() + "\t" + line_sets[6].strip() + "\t" + line_sets[7].strip()
                            sql = sql + "\t" + line_sets[8].strip()
                            sql = sql + "\t" + line_sets[int(patientLOC)].strip() + "\t" + line_sets[int(maternalLOC)].strip()
                            fileName = outDir + "/out_" + fid + ".vcf"
                            relativeLOC = extractRelative(fid)
                            if len(relativeLOC) != 0:
                                sql1 = line_sets[int(relativeLOC[0] + 8)].strip()
                                for i in range(1, len(relativeLOC)):
                                    sql1 = sql1 + "\t" + line_sets[int(relativeLOC[i] + 8)].strip()
                                sql = sql + "\t" + sql1
                            with open(fileName, "a") as f:
                                f.writelines(sql + "\n")
                    else:
                        pass
                else:
                    #print(line)
                    if int(patientDP) >= 10 and int(maternalDP) >= 10 and int(paternalDP) >= 10:
                        if (patientGT == "0/1" or patientGT == "0|1" or patientGT == "1/1" or patientGT == "1|1" or patientGT == "1/0" or patientGT == "1|0") and (maternalGT == "0/0" or maternalGT == "0|0") and (paternalGT == "0/0" or paternalGT == "0|0"):
                            sql = line_sets[0].strip() + "\t" + line_sets[1].strip() + "\t" + line_sets[2].strip() + "\t" + line_sets[3].strip()
                            sql = sql + "\t" + line_sets[4].strip() + "\t" + line_sets[5].strip() + "\t" + line_sets[6].strip() + "\t" + line_sets[7].strip()
                            sql = sql + "\t" + line_sets[8].strip()
                            sql = sql + "\t" + line_sets[int(patientLOC)].strip() + "\t" + line_sets[int(maternalLOC)].strip() + "\t" + line_sets[int(paternalLOC)].strip()
                            fileName = outDir + "/out_" + fid + ".vcf"
                            relativeLOC = extractRelative(fid)
                            if len(relativeLOC) != 0:
                                sql1 = line_sets[int(relativeLOC[0]+8)].strip()
                                for i in range(1,len(relativeLOC)):
                                    sql1 = sql1 + "\t" + line_sets[int(relativeLOC[i]+8)].strip()
                                sql = sql + "\t" + sql1
                            with open(fileName,"a") as f:
                                f.writelines(sql+"\n")
                    else:
                        pass
