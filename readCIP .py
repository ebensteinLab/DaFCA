#!/usr/bin/env python
# Code Written by Jeffrey G. Reifenberger for BioNano Genomics
import sys
import re
import math
import argparse
import os.path
import numpy as np
import random
import struct as st

def checkFile(bnxFile,cipxFile):

    if len(bnxFile) == 0:
        print('##### No .bnx file given! (-b)')
        sys.exit()

    if len(cipxFile) == 0:
        print('##### No .cipx file given! (-x)')
        sys.exit()




#this code is for debugging purposes. Commented in normal code.
def compareLength(originalToMolID,molPosHash,bpp):
    #outfile name
    ofile = open('bnxVScipx_lengthDiff.txt','w')

    #tempFOV = [fov#,xStart,xEnd,yStart,cipFilePos,lenFOVPx,ystitchOffSet]

    ofile.write('oMolID\tbnxLength\tcipxLength\tDiff\n')

    #list of original moleculeIDs.
    okeys = list(originalToMolID.keys())

    #iterate over all moleculeIDs.
    for o in okeys:
        #write orig. mol ID, and length of molecule from .bnx file.
        ofile.write(str(o)+'\t'+str(originalToMolID[o][1]))

        lenPix = 0
        for i in range(0,len(molPosHash[o])):
            lenPix+=molPosHash[o][i][5] - molPosHash[o][i][6]
            #compute length of molecule from .cipx file in pixels.
            #molPosHash[o][i][5] -> length of molecule segment in FOV
            #molPosHash[o][i][6] -> overlap length between adjacent FOVs

        #convert lenPx to bp, compute difference between .bnx file and .cipx file.
        ofile.write('\t'+str(lenPix*bpp)+'\t'+str(originalToMolID[o][1] - lenPix*bpp)+'\n')

    ofile.close()

    #exit program.
    sys.exit()

        



def readBNXFile(filename):
    #Functions main purpose is to convert from original moleculeID in the .cipx file to the moleculeID in the .bnx file.
    #original molecule ID: 7th entry in 0 line, moleculeID: 2nd entry in 0 line of .bnx file
    #all alignments from refaligner are based on moleculeID, not original moleculeID.
    #need to read all continuous profiles and replace with proper moleculeID
    infile = open(filename,'rU')

    #originalToMolID: hash, key = originnal moleculeID, value = [moleculeID, moleculelength(bp), runID, and column]
    originalToMolID = {}

    for line in infile:
        if line[0] == '#':  #skip # header lines
            continue
        elif line[0] == '0':    #only read 0 line from each molecule (ignore 1, 2, and QX codes) 
            word = line.split()

            if str(word[6]) in originalToMolID:
                #double check that original molecule ID is unique in .bnx file. Should be unique, but just in case exit code...
                print(str(word[6])+' already exists in bnx file...')
                sys.exit()
            else:
                originalToMolID[str(word[6])] = [str(word[1]),float(word[2]),int(word[11]),int(word[12])]
                        #word[6] = originalMoleculeID, word[1]=molID,word[2]=molLength,word[11]=runID,word[12]=column

    infile.close()
    return originalToMolID

def readCIPXFile(pColor,filename,prefix):
    infile = open(filename,mode='rU')

    #molPosHash, hash. key->original molID, value->list of lists for each segment of molecule
    #within each sublist is a fov #, x, y, coordinates of segment, and start location of profile in
    #.cip file.
    molPosHash = {} 

    #cipNameHash, hash. key->.cip filename, value is list of lists->[[original molecule ID_1, bank ID (1,2,3, or 4)], [original molID_2, bankID]... etc]
    #list of list for each is all the original molecule IDs in .cip file.
    cipNameHash = {}

    for line in infile: 
        if line[0] == '#' or line[0] == 'M' or line[0] == 'O':
            #skip headers in .cipx file.
            continue
        elif len(line) == 1:
            #skip blank line
            continue
        else:
            word = line.split()

            #split file name provided in .cipx to obtain bankID
            temp = word[1].split('_')   
            bankID = temp[1]

            #generate .cip name for color given by user (pColor).
            cipName = str(word[1])+'_'+str(pColor)+'.cip'

            #store .cip filename of interest with original moleculeID and bankID in cipNameHash
            if cipName in cipNameHash:
                cipNameHash[cipName].append([str(word[0]),bankID])  #store every molecule for each .cip file as list of lists
            else:
                cipNameHash[cipName] = [[str(word[0]),bankID]]

            #number of FOVs 
            numFOVs = int(word[2])

            #starting FOV of molecule
            fovStart = int(word[3])

            #Compute the number of FOV entries for the molecule.
            numEntry = float(len(word) - 4)/6.0

            #tempFOV, store each molecule segment info before passing to tempFull
            tempFOV = []
            #tempFull, store each molecule segments info before passing to molPosHash
            tempFull = []

            #Fov info starts at w = 4 entry in .cipx file (zero based)
            w = 4
            for i in range(0,numFOVs):  #interate through each FOV of molecule.
                tempFOV.append(fovStart)   #starting fov number placed in tempFOV. 
                for j in range(w,w+6):  #interate through individual molecule segments.
                    tempFOV.append(float(word[j]))

                w += 6  #add 6 to w to move to next segment.

                #add tempFOV to tempFull for molecule segment.
                #tempFOV = [fov#,xStart,xEnd,yStart,cipFilePos,lenFOVPx,ystitchOffSet]
                #yEnd (not in .cipx file) can be computed from xStart, xEnd, yStart, and lenFOVPx
                #location of start of profile trace in .cip file is cipFilePos
                #ystitchOffset allows user to know overlap region (in pixels) between two adjacent FOVs
                tempFull.append(tempFOV)
                
                #reset tempFOV for next segment
                tempFOV = []
                fovStart += 1   #iterate to next FOV #.

            molPosHash[str(word[0])] = tempFull #store all segments of molecule in molPosHash
            tempFull = []   #reset tempFull for next molecule.

    infile.close()
    return molPosHash,cipNameHash


def readCIPFile(outfile,filename,cipNameHash,molPosHash,originalToMolID,bpp):

    #open .cip file.
    with open(filename, mode='rb') as file:
        filecontent = file.read()

    #iterate through each molecule, c, is the .cip file
    for c in cipNameHash[filename]:
        m = c[0]       #m is original molecule ID
        bankID = c[1]   #bankID for molecule
        
        #tempFOV stores each intensity value for the molecule segment
        tempFOV = []
        #tempFull is list of lists, stores each segment's intensity values. Eventually printed to _profileTrace.txt file.
        tempFull = []

        #iterate through each segment list in molPosHash[m]
        #Each segment list in molPosHash contains: [fov#,xStart,xEnd,yStart,cipFilePos,lenFOVPx,ystitchOffSet]
        for i in range(0,len(molPosHash[m])):
            cipStart = int(molPosHash[m][i][4])     #location of start byte of profile trace in .cip file.

            ####IMPORTANT#####
            #the profile traces for each segment of a molecule in the .cip file should be flanked by zero. This allows
            #the user to double check if the profile grabed from .cip file is correct.
            #cipStart - 4 should be 0
            #cipStart + 4 *(lenFOVPx) should be 0 (remember that python indexes at 0, not 1)
            for j in range(0,int(molPosHash[m][i][5])):   #iterate through the length of the segment
                tempFOV.append(round(st.unpack('<f',filecontent[cipStart:cipStart+4])[0],3))    #round intensity values to 3 decimal places
                #st.unpack returns a tuple in which the 2nd value is empty, e.g (#,). Hence why append [0] value of tuple to tempFOV
                #one pixel of the profile is equivalent to 4 bytes in .cip file, values are little endian and should be converted to 32 bit float.
                cipStart += 4   #iterate cipStart 4 positions to read the next camera pixel intensity value.

            tempFull.append(tempFOV)    #add tempFOV to tempFull
            tempFOV = []    #reset tempFOV
        
        #print profile trace in .txt file.
        printProfileTrace(outfile,m,bankID,originalToMolID,molPosHash,tempFull,bpp)
        tempFull = []


def computeYend(x_s,x_e,y_s,L):
    #compute y_end from quadratic equation.
    return round(y_s+np.sqrt(L**2 - (x_e-x_s)**2),1)

def computeMolStitch(molLen,molList,bpp):
    #assume bpp is given by user. bpp from .bnx file. typically bpp = 375.0
    #store full x,y values for each segment.
    xyList = []

    #molist contains molecule segment positions from molPosHash -> [fov#,xStart,xEnd,yStart,cipFilePos,lenFOVPx,ystitchOffSet]

    #iterate through each segment of molecule. calculate y_end value for segment.
    for i in range(0,len(molList)):
        y_end = computeYend(molList[i][1],molList[i][2],molList[i][3],molList[i][5])

        #add y_end to list of x,y coordinates.
        #[fov#,xStart,yStart,xEnd,yEnd]
        xyList.append([molList[i][0],molList[i][1],molList[i][3],molList[i][2],y_end])

    #the bpStart position for the first molecule segment is at 0bp
    bpStartFOV = [0]
    lenPx = molList[0][5]   #length of segment in pixels of first segment
    #lenPx is accumlative length of molecule after each molecule segment. Takes into account overlaping stitch regions between FOV.
    for i in range(1,len(molList)): #start at 2nd segment of molecule.
        bpStartFOV.append(bpp*(lenPx - molList[i][6]))  #start position of next segment in bp on molecule.
        #subtrace yStitchOffSet value from accumlated length of molecule, lenPx

        #new lenPx value is equal to segment length of ith segment of molecule minus overlap region
        lenPx += molList[i][5] - molList[i][6]  
        
    return bpStartFOV,xyList
        

    
def printProfileTrace(outfile,o,bankID,originalToMolID,molPosHash,tempFull,bpp):
    #outfile is object for _profileTrace.txt file.
    #o is original molecule ID
    #bankID is bank of molecule o
    #originalToMolID, hash. key originalMolID o, value is list [molID, molLength,runID,column]
    #molPosHash, hash. key originalMolID o, value is list of lists. Each sublist contains coordinates for segment.
    #sublist in molPosHash = [fov#,xStart,xEnd,yStart,cipFilePos,lenFOVPx,ystitchOffSet]
    #tempFull: List of lists, each sublist contains intensity values per pixel of molecule.

    #Start each molecule info line with 'm'
    #write moleculeID, molecule Length (bp), runID, column, bankID
    outfile.write('m\t'+str(originalToMolID[o][0])+'\t'+str(originalToMolID[o][1])+'\t'+str(originalToMolID[o][2])+'\t'+str(originalToMolID[o][3])+'\t'+str(bankID))

    #compute bpFOV start position for each segment in molecule. xyList contains full x,y coodinates for each segment.
    bpStartFOV,xyList = computeMolStitch(originalToMolID[o][1],molPosHash[o],bpp)

    #iterate over xyList values. write them to 'm' line of _profileTrace.txt
    for i in range(0,len(xyList)):
        for j in range(0,len(xyList[i])):
            outfile.write('\t'+str(xyList[i][j]))
    outfile.write('\n') #new line

    #write #fov number of segment, bpStart location in molecule of respective segment.
    for i in range(0,len(tempFull)):
        outfile.write(str(molPosHash[o][i][0])+'\t'+str(bpStartFOV[i]))

        #write itensity value of each pixel for molecule segment
        for j in range(0,len(tempFull[i])):
            outfile.write('\t'+str(tempFull[i][j]))
        outfile.write('\n') #new line.

        #repeat for each molecule segment

    
        


######################################################################

parser = argparse.ArgumentParser(description='Reads .cip file in order to obtain the continuous profile values for each molecule. Data is output to a text file.')
parser.add_argument("-b", "--bnxFile", help=".bnx filename",type=str,default='')
parser.add_argument("-i", "--profileColor", help="channel number (color) of .cip file. Values can be 1, 2, or 3. 1 is always blue. 2 is typically red. 3 is typically green. Double check to be safe. default=3",type=int,default=3)
parser.add_argument("-x", "--cipXFile", help=".cipx filename",type=str,default='')
parser.add_argument("-c", "--bppValue", help="bpp value of respective .bnx file, should be 375. default=375",type=float,default=375.0)
parser.add_argument("-p", "--prefix", help="prefix for output data files",type=str,default='some_great_data')
args = parser.parse_args()

#Check that user gives bnx and cipx files.
checkFile(args.bnxFile,args.cipXFile)

print('#### ProfileColor Number (-i):\t'+str(args.profileColor)+'\n')

print('Reading .bnx file\t'+args.bnxFile)

# moleculeID values in the .cipx file are the original moleculeID from the .bnx file. Not the moleculeID
originalToMolID = readBNXFile(args.bnxFile)

print('Reading .cipx file\t'+args.cipXFile)
#read .cipx file. This contains x,y pixel values for each segment of the molecule per FOV, and position in .cip file
#to find continuous profile in color of interest.
molPosHash,cipNameHash = readCIPXFile(args.profileColor,args.cipXFile,args.prefix)

#function compareLength compares the total length of each molecule in the .bnx file to that calculated from the .cipx file.
#this function is commented in a regular run, but potential useful for debugging.
#compareLength(originalToMolID,molPosHash,args.bppValue)


#ckeys is list of all .cip filenames of profile color of interest.
ckeys = list(cipNameHash.keys())
ckeys.sort()        #sort by filename.


#Header in output file prefix_profileTrace.txt
outfile = open(args.prefix+'_profileTrace.txt','w')
outfile.write('#m\tmolID\tmolLength(bp)\trunID\tcolumn\tBankID\tfov[1]\txStart[1]\tyStart[1]\txEnd[1]\tyEnd[1]')
outfile.write('\tfov[2]\txStart[2]\tyStart[2]\txEnd[2]\tyEnd[2]')
outfile.write('\tfov[3]\txStart[3]\tyStart[3]\txEnd[3]\tyEnd[3]')
outfile.write('\tfov[4]\txStart[4]\tyStart[4]\txEnd[4]\tyEnd[4]\n')    
outfile.write('#fovNum\tbpStart\tprofileTraceIntensity -->\n')

#iterate through each .cip file reading all the profile traces from each molecule. Plase output in outfile.
for c in ckeys:
    print('Reading .cip file:\t'+c)
    readCIPFile(outfile,c,cipNameHash,molPosHash,originalToMolID,args.bppValue)


outfile.close() 




