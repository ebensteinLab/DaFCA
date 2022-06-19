#!/usr/bin/env python
# Code Written by Jeffrey G. Reifenberger for BioNano Genomics
import sys
import re
import math
import argparse
import os.path
import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def checkFile(bFile,pFile,xFile,cFile):

    if len(bFile) == 0:
        print('#### no bnx file given! (-b)')
        sys.exit()

    if len(pFile) == 0:
        print('#### no profile data file given! (-m)')
        sys.exit()

    if len(xFile) == 0:
        print('#### no xmap file given! (-x)')
        sys.exit()        

    if len(cFile) == 0:
        print('#### no cmap file given! (-x)')
        sys.exit()        

#read region string provided by user. Split by comma.
def parseRegion(rString):
    return rString.split(',')

#read cmap file provided by user.
def readCmapFile(filename):
    infile = open(filename,'rU')

    #cmapHash is a hash that stores the nick/DLS site in bp.  The key for cmapHash is the chrm,label#
    cmapHash = {}

    for line in infile:
        if line[0] == '#':
            continue
        else:
            word = line.split()

            if int(word[4]) == 1:   #if word[4] == 0, reached end of the chromosome. if word[4] = 0, then word[5] = length of chromosome.
                cmapHash[str(word[0])+','+str(word[3])] = float(word[5])
                #word[0] = chromosome number, word[3] = label number (nick/DLS), word[5] = site of nick/DLS in bp

    infile.close()
    return cmapHash


#read xmapfile, selecting for alignments within regionList.
def readXmapFile(regionList,filename):
    infile = open(filename,'rU')

    #xMapHash, key=moleculeID, values=[chromosomeID, direction of alignment(+/-), alignment -logPvalue (confidence score)]
    xMapHash = {}

    #alignStrHash, key=moleculeID, value = [[refPos1, mol_label1],[refPos2,mol_label2],...[refPosN,mol_labelN]] => list of lists.  
    #Each sub-list contains the reference position from the cmap file that a particular label aligned to from the molecule.
    alignStrHash = {}

    for line in infile:
        if line[0] == '#':
            continue
        else:
            word = line.split()
            # word[2] = chromosome number from molecule in xmap file, regionList[0] is chromosome number given by user.
            if int(word[2]) == int(regionList[0]):
                temp = re.findall(r'(\d+),(\d+)',str(word[13])) #split word[13], alignment string, into list of lists, store in alignStrHash

                #check that alignment region within user values given in regionList
                if int(temp[0][0]) >= int(regionList[1]) and int(temp[-1][0]) <= int(regionList[2]):
                    xMapHash[str(word[1])] = [int(word[2]),str(word[7]),float(word[8])]
                                #word[1] = molID, word[2]=chrmID,word[7]=dir,word[8]=confScore (-logPvalue)

                    #temp = [(refPos1,labelPos1), (refPos2,labelPos2),....(refPosN,labelPosN),]
                    alignStrHash[str(word[1])] = temp

    infile.close()
    return xMapHash,alignStrHash

#read profile trace and store data in hashes.
def readProfileFile(xMapHash,filename):
    infile = open(filename,'rU')

    #molStats hash: key = molID, value = list of the general (e.g. length, x pos, y pos, etc.) stats of the molecule.
    molStats = {}
    #bpStart hash: key = molID, value = list of start values in bp for each FOV of the profile trace.
    bpStart = {}
    #profileTrace hash: key = molID, value = list of lists. Each sublist is a profile trace from an FOV.
    profileTrace = {}


    for line in infile:
        if line[0] == '#':
            continue
        elif line[0] == 'm':
            word = line.split()
            keep = 0

            molID = str(word[1])
            if str(word[1]) in xMapHash:    #check if moleculeID, word[1] in xMapHash. if yes, store values in appropriate hash.
                molStats[str(word[1])] = word[2:]
                # word[1] = molID
                #[word[2]=molLength,word[3]=runID,word[4]=column,word[5]=fovStart,word[6]=xStart,word[7]=yStart,word[8]=fovEnd,word[9]=xEnd,word[10]=yEnd]
                keep = 1
            else:
                keep = 0

        else:
            if keep == 1:
                word = line.split()

                #store bpStart values in hash
                if molID in bpStart:
                    bpStart[molID].append(float(word[1]))
                else:
                    bpStart[molID] = [float(word[1])]

                #store profile trace values in hash.
                if molID in profileTrace:
                    profileTrace[molID].append([float(x) for x in word[2:]])    #converts all data in the profile trace to float
                else:
                    profileTrace[molID] = []
                    profileTrace[molID].append([float(x) for x in word[2:]])

    infile.close()
    return molStats,bpStart,profileTrace

#reads bnx file, returns the label position and intensity per molecule
def readBNXfile(xMapHash,filename):
    infile = open(filename,'rU')

    #bnxData hash: keys= molID,1 or molID,QX12. Values are the label positions in bp (1) or the label intensity values (QX12)
    bnxData = {}

    for line in infile:
        if line[0] == '#':
            continue
        elif line[0] == '0':
            word = line.split()
            keep = 0

            #keep molecule if moleculeID=word[1] is in xMapHash, which is the region of interest given by user.
            if str(word[1]) in xMapHash:
                molID = str(word[1])
                keep = 1
            else:
                keep = 0

        #if keep = 1, then keep the values of label locations in bp.
        elif line[0] == '1' and keep == 1:
            word = line.split()
            tempL = [float(x) for x in word[1:]]

        #if keep = 1, then keep the value of the label intensity for each site.
        elif line[0:4] == 'QX12' and keep == 1:
            word = line.split()
            tempSNR = [float(x) for x in word[1:]]

            #store the data in bnxData.
            bnxData[molID+',1'] = tempL
            bnxData[molID+',QX12'] = tempSNR

    infile.close()
    return bnxData

#assign every label on a molecule to a location on the reference based on alignment
def organizeBNXData(bnxData,xMapHash,alignStrHash,cmapHash):
    mkeys = list(alignStrHash.keys())     #mkeys is list of molecule IDs that align to the region of interest given by user.

    #molAlignHash hash: key = moleculeID, value = list of lists. each sublist is per label [[label# molecule, bp value of label in molecule, intensity of label,
    #location of alignment to reference in bp from _r.cmap file ], .... continue for each label ...]
    molAlignHash = {}

    for m in mkeys:
        chrmID = str(xMapHash[m][0])
        temp = []   #temp temporarily stores list of lists before passing into molAlignHash.
        for i in range(0,len(alignStrHash[m])):
            refN = str(alignStrHash[m][i][0])       # reference siteID from cmap file
            labelN = int(alignStrHash[m][i][1]) - 1 # label number from molecule, covert to zero base for python. 
            k = chrmID+','+refN                     #k = key for cmapHash, chromosomeID,nick/DLS site ID
            temp.append([labelN,bnxData[m+',1'][labelN],bnxData[m+',QX12'][labelN],cmapHash[k]])
                        #[label # on molecule, bp value of label on molecule, intensity value of label on molecule, location in bp of alignment to cmap]

        temp.append([-1,bnxData[m+',1'][-1],-1,-1])  #last entry is the length of molecule, other values in list set to -1.
        molAlignHash[m] = temp      #past the temp list to molAlignHash

    return molAlignHash

#Compute the starting point of the profile in the reference based on label alignment positions in the reference.
def findStartPoint(m,fovStart,startOffset,refLabelStart):
    #m is molecule ID
    #fovStart: list of starting points in bp of the molecule for each FOV
    #startOffset: This value is based on the fact that labels that align on the molecule are often not a the left/right end of a molecule but rather in their interior. 
    #account for this overhang offset from label to end of molecule.
    #refLabelStart, alignment location in bp, of the first label in a molecule to the reference.

    if fovStart == 0.0:
        return refLabelStart - startOffset
    else:
        return refLabelStart - startOffset + fovStart


def addXBP(measBPP,startFOV,mTraceList):

    #temp: list of lists. Each sublist contains the bp of the profile based on the alignment to the reference and the profile value
    #[[refBP_1,profileValue_1], [refBP_2,profileValue_2], [refBP_3,profileValue_3] .... [refBP_N,profileValue_N]]
    temp = []

    for j in range(0,len(mTraceList)):
        temp.append([startFOV+j*measBPP,mTraceList[j]])

    return temp

#algorithm reverses fov start sites for molecule that align in the - direction.
def computeReverseBPStart(fovStartList,molLength,tempTraceList):
    #list of the new start points, in bp, for each FOV if the molecule aligns in the reverse, -, direction
    rev_bpStart = []

    for j in range(0,len(fovStartList)):
        rev_bpStart.append(round( molLength - (fovStartList[j]+375.0 * len(tempTraceList[j])),1 ) )
        #the new start site for each FOV for reverse alignment = molecule_length - (fovStart_j + 375 * lengthFOV_px_j)
        #e.g. molecule length = 319,875bp with start sites at 0 and 92,625bp and a length of FOV_0 = 454px  and FOV_1 = 606px
        #Then the new start points will be 319,875bp - (0 + 375*454) = 149,625bp and 319,875bp - (92,625 + 375*606) = 0

    rev_bpStart.sort()  #values are sorted to go from lowest to highest, the conversion results in the opposite

    return rev_bpStart

def computeReverseTrace(tempTrace):

    #rev_Trace: list of lists, place reversed profile traces
    rev_Trace = []

    for j in range(len(tempTrace)-1,-1,-1): #reverse order of profile trace.
        rev_Trace.append(tempTrace[j][::-1])    #flip profile trace for each fov

    return rev_Trace


def organizeProfileTraceData(measBPP,molAlignHash,xMapHash,bpStart,profileTrace):
    mkeys = list(molAlignHash.keys())     #mkeys, list of molecule IDs that align to the region of interest.

    #profileAlignHash: hash, key = moleculeID, value = list of lists [[bpReference_1,profileValue_1], [bpReference_2,profileValue_2],... [bpReference_N,profileValue_N]]
    profileAlignHash = {}

    for m in mkeys:
        #if the alignment is reversed, then do as follows
        if xMapHash[m][1] == '-':
            #compute the bp start sites for each FOV based on reverse alignment.
            rev_bpStart = computeReverseBPStart(bpStart[m],molAlignHash[m][-1][1],profileTrace[m])
                            #bpStart = values from morse.Profile, molAlignHash[m][-1][1] = molecule Length, profileTrace[m] = list of lists, profile for each fov.

            #flip the profile traces for each FOV.
            rev_Trace = computeReverseTrace(profileTrace[m])

            for i in range(0,len(rev_bpStart)):
                #compute the start point for each FOV in the reference .cmap bp frame
                #convert from 375 bp/px from the .bnx file to the actual measured bp/px from alignment, measBPP
                startFOV_ref = findStartPoint(m,(measBPP/375.0)*rev_bpStart[i],(measBPP/375.0)*(molAlignHash[m][-1][1] - molAlignHash[m][0][1]),molAlignHash[m][0][3])
                #molAlignHash = [label # on molecule, bp value of label on molecule, snr value of label on molecule, location in bp of alignment to cmap]

                #after new starting points have been calculated and profile traces have been flipped, add the bp values in the reference frame.
                tempTrace = addXBP(measBPP,startFOV_ref,rev_Trace[i])

                #add tempTrace to profileAlignHash
                if m in profileAlignHash:
                    profileAlignHash[m].append(tempTrace)
                else:
                    profileAlignHash[m] = []
                    profileAlignHash[m].append(tempTrace)
        #if the alignment is not reversed, then go here. No need to reverse any of the profile values or start locations per FOV.
        else:
            for i in range(0,len(bpStart[m])):
                #See the comments for reverse aligment above. Same thing, but this time no reverse of profile values or recalculating of FOV start sites.
                startFOV_ref = findStartPoint(m,(measBPP/375.0) * bpStart[m][i],(measBPP/375.0)*molAlignHash[m][0][1],molAlignHash[m][0][3])

                tempTrace = addXBP(measBPP,startFOV_ref,profileTrace[m][i])

                if m in profileAlignHash:
                    profileAlignHash[m].append(tempTrace)
                else:
                    profileAlignHash[m] = []
                    profileAlignHash[m].append(tempTrace)

    return profileAlignHash

#generate list of label positions and intensity values for plotting.
def generateLabelLists(labelAlignList,offSetL):
    x = []
    y = []

    for j in range(0,len(labelAlignList)-1):    #last value is molecule length
        x.append(round(labelAlignList[j][3]/1.0e6,3))
        y.append(3000.0+offSetL)    #intensity values are arbitrary. merely for plotting on same graph.

    return x,y

#generate list of label positions and intensity values for data text file.
def generateLabelsForFile(labelAlignList):
    x = []
    y = []

    for j in range(0,len(labelAlignList)-1):    #last value is molecule length, hence why -1.
        x.append(round(labelAlignList[j][3]/1.0e6,3))
        y.append(labelAlignList[j][2])  #intensity value of label.

    return x,y,labelAlignList[-1][1]

#generate lists of profile values for plotting.
def generateProfileList(profileAlignList,offSetM):
    x = []
    y = []

    for j in range(0,len(profileAlignList)):
        x.append(round(profileAlignList[j][0]/1,6))
	#x.append(profileAlignList[j][0])
        y.append(profileAlignList[j][1]+offSetM)    #profile value is intensity with offset so traces can be stacked.

    return x,y

#generate lists of profile values for data text file.
def generageProfileforFile(profileX_file,profileY_file,profileAlignList):
    x = []
    y = []

    for j in range(0,len(profileAlignList)):
        x.append(round(profileAlignList[j][0]/1,6))
	#x.append(profileAlignList[j][0])
        y.append(profileAlignList[j][1])    #profile value is intensity, no offset.

    profileX_file.append(x)
    profileY_file.append(y)

    return profileX_file,profileY_file


def plotData(regionList,molAlignHash,profileAlignHash,xMapHash,prefix):

    #colorList: list of profile color options.
    colorList = ['#0316FB','#2FAEB9','#7C1FD9','#6ED91F']
    #generate pdf file with pages
    outPDF = PdfPages(prefix+'_labelAlign_ProfileTrace.pdf')

    #mkeys = list of moleculeID keys.
    mkeys = list(molAlignHash.keys())

    print('Number of molecules that align to ROI:',len(mkeys))

    #offSetL = offset of label y values.
    offSetL = 40000
    #offSetM = offset of profile y values.
    offSetM = 0
    #open figure, label axis
    plt.figure()
    plt.title('Profile Traces with Labels to Chrm '+str(regionList[0])+' refSite# '+str(regionList[1])+' to '+str(regionList[2]))
    plt.xlabel('Reference Alignment (Mbp)')
    plt.ylabel('Intensity (counts)')
    legendList = ['Label','FOV0','FOV1']

    #open a .txt file that the raw data per molecule is written too. this will allow user to plot data per molecule if desired.
    outfile = open(prefix+'_labelAndProfileRawData.txt','w')
    outfile.write('#molID\tmolLength(bp)\tnumLabels\n')
    outfile.write('#label_bp\tvalues--->\n')
    outfile.write('#label_intensity\tvalues--->\n')
    outfile.write('#profile_bp_0\tvalues--->\n')
    outfile.write('#profile_intensity_0\tvalues--->\n')
    outfile.write('#profile_bp_N\tvalues--->\n')
    outfile.write('#profile_intensity_N\tvalues--->\n')


    #interate through each aligned molecule
    for m in mkeys:
        #generate x and y value lists of label positions for plotting.
        labelX,labelY = generateLabelLists(molAlignHash[m],offSetL)

        #generate x and y value lists of label positions for text data file.
        labelX_file,labelY_file,molLength = generateLabelsForFile(molAlignHash[m])

        #plot label positions .
        plt.plot(labelX,labelY,marker='o',color='#D9951F',linestyle='none',markersize=1)

        profileX_file = []
        profileY_file = []
        #plot profile trace for each FOV in the molecule
        for i in range(0,len(profileAlignHash[m])):
            #generate x and y value lists for profile data.
            profileX,profileY = generateProfileList(profileAlignHash[m][i],offSetM)

            #generate profile values for data text file.
            profileX_file,profileY_file = generageProfileforFile(profileX_file,profileY_file,profileAlignHash[m][i])

            #plot profile data as a line.
            plt.plot(profileX,profileY,color=colorList[i],linewidth=1)

        #interate offSet values for each molecule so they can stack in the figure.
        offSetL += 500
        offSetM += 2000

        #print molecule data in a .txt file as well for user to plot per molecule.
        printMoleculeData(m,molLength,outfile,labelX_file,labelY_file,profileX_file,profileY_file)

    plt.legend(legendList)
    outPDF.savefig(dpi = 1000)
    plt.close()

    outPDF.close()
    outfile.close()

def printMoleculeData(m,molLength,outfile,labelX_file,labelY_file,profileX_file,profileY_file):
    #outfile has a series of rows for each molecule.
    #each row for molecule begins with a text identifer.
    #molID: contains molID from .bnx file, molecule length, number of labels
    #label_bp: label bp values from alignment to reference
    #label_intensity: label intensity values.
    #profile_bp_fov#: profile bp value from alignment to reference
    #profile_intensity_fov#: profile intenisty values.
    #repeat for each molecule that aligns to location given by user.

    outfile.write('molID\t'+str(m)+'\t'+str(molLength)+'\t'+str(len(labelX_file))+'\n')

    outfile.write('label_bp')
    for x in labelX_file:
        outfile.write('\t'+str(x))
    outfile.write('\n')
        
    outfile.write('label_intensity')
    for y in labelY_file:
        outfile.write('\t'+str(y))
    outfile.write('\n')

    for i in range(0,len(profileX_file)):
        outfile.write('profile_bp_'+str(i))
        for x in profileX_file[i]:
            outfile.write('\t'+str(x))
        outfile.write('\n')

        outfile.write('profile_intensity_'+str(i))
        for y in profileY_file[i]:
            outfile.write('\t'+str(y))
        outfile.write('\n')        
            


parser = argparse.ArgumentParser(description='Sorts through alignment of label data and correlates to profile data. Plots coverage for a region of interest based on user input.')
parser.add_argument("-b", "--bnxFile", help="name of bnx file",type=str,default='')
parser.add_argument("-m", "--profileFile", help="name of profile trace file.",type=str,default='')
parser.add_argument("-x", "--xmapFile", help="name of .xmap file.",type=str,default='')
parser.add_argument("-c", "--cmapFile", help="name of _r.cmap file.",type=str,default='')
parser.add_argument("-r", "--region", help="region of interest. format chrmID,nickStart,nickEnd",type=str,default='1,100,125')
parser.add_argument("-a", "--mBPP", help="measured bpp value. default=375",type=float,default=375.0)
parser.add_argument("-p", "--prefix", help="prefix for output data files",type=str,default='some_great_data')
args = parser.parse_args()

#check user gave correct files.
checkFile(args.bnxFile,args.profileFile,args.xmapFile,args.cmapFile)

#read and parse the region of interest string input by the user.
#regionList, list contains chromosome number, start site, end site.
regionList = parseRegion(args.region)

#Read _r.cmap file. 
print('Reading cmap File:\t',args.cmapFile)
cmapHash = readCmapFile(args.cmapFile)

#read .xmap file
print('Reading xmap File:\t',args.xmapFile)
xMapHash,alignStrHash = readXmapFile(regionList,args.xmapFile)

#read profile file, select molecules that are in xMapHash, that align to region of interest.
print('Reading profile trace File:\t',args.profileFile)
molStats,bpStart,profileTrace = readProfileFile(xMapHash,args.profileFile)

#read bnx file. select molecules in xMapHash that are in the region of interest.
print('Reading bnx File:\t',args.bnxFile)
bnxData = readBNXfile(xMapHash,args.bnxFile)

#assign each label in the molecule with a reference position.
print('Matching Labels to Reference:')
molAlignHash = organizeBNXData(bnxData,xMapHash,alignStrHash,cmapHash)

#Align the profile data to the reference based on the alignment of the labels.
print('Matching Profile Traces with Alignment:')
profileAlignHash = organizeProfileTraceData(args.mBPP,molAlignHash,xMapHash,bpStart,profileTrace)

#generate a simple coverage plot for the labels and the profile data for the region of interest.
plotData(regionList,molAlignHash,profileAlignHash,xMapHash,args.prefix)



