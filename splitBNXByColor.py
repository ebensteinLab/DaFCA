#!/usr/bin/env python
# Code Written by Jeffrey G. Reifenberger for BioNano Genomics
import sys
import re
import math
import argparse
import os.path
import numpy as np
import random


#Function findColorCode determines the color of the nickase line in the .bnx file. 
#There are a series of potential color names in the .bnx file that the function searches.
def findColorCode(line,colorHash):
    #split # Nickase Recognition Site by white space.
    word = line.split()

    #color number of line is the 5th value (0 start) of the list. 0:1 is number, 1:2 is :
    colorNum = str(int(word[4][0:1]))   #colorNum for the labels will be either 1 or 2.

    #nick sequence/color is last entry in word.
    # split squence/color by split with ;
    nickColor = word[-1].split(';')
    #find color by splitting last value (1) of nickColor list with '_'. Sometimes color is Green_001
    color = nickColor[1].split('_')
    #color name is first entry of color.
    colorName = color[0]

    #search through potential .bnx color options. List is not extensive, but should be largely correct.
    #place values in colorHash: key = color number (1 or 2), value is either Red or Green depending on search.
    if colorName == 'bngflrd001':
        colorHash[colorNum] = 'Red'
    elif colorName == 'BNGFLRD001':
        colorHash[colorNum] = 'Red'
    elif colorName == 'red':
        colorHash[colorNum] = 'Red'
    elif colorName == 'bngflgr001':
        colorHash[colorNum] = 'Green'
    elif colorName == 'BNGFLGR001':
        colorHash[colorNum] = 'Green'
    elif colorName == 'green':
        colorHash[colorNum] = 'Green'
    else:
        colorHash[colorNum] = colorName

    return colorNum,colorHash

#splitBNXFile reads input bnx file and splits into green or red bnx file depending on color scheme.
def splitBNXFile(filename,prefix):
    #open input .bnx file
    infile = open(filename,'r')

    #open a green and red .bnx file for output.
    gfile = open(prefix+'_GreenOnly.bnx','w')
    rfile = open(prefix+'_RedOnly.bnx','w')


    #colorHash: stores color information from bnx file. key is 1 or 2. value is green or red.
    colorHash = {}

    for line in infile:
        #for hashed lines in input .bnx file, they need to be slightly altered for a single color .bnx file.
        if line[0] == '#':
            #original .bnx file has 2 color (green and red), the new .bnx files will have a single color (green or red).
            if line[0:17] == '# Label Channels:':
                rfile.write('# Label Channels:\t1\n')
                gfile.write('# Label Channels:\t1\n')

            #colorHash stores color information from bnx file. key is 1 or 2, value Green or Red
            elif line[0:26] == '# Nickase Recognition Site':
                colorNum,colorHash = findColorCode(line,colorHash)

                #Depending on results from findColorCode, write new lines in split .bnx files as either green or red.
                if colorHash[str(colorNum)] == 'Green':
                    word = line.split('\t')
                    gfile.write('# Nickase Recognition Site 1:\t'+word[1])
                elif colorHash[str(colorNum)] == 'Red':
                    word = line.split('\t')
                    rfile.write('# Nickase Recognition Site 1:\t'+word[1])

            #split .bnx files are now single color. Write down the #1 lines.
            elif line[0:2] == '#1':
                gfile.write(line)
                rfile.write(line)
            #split .bnx file are now single color, ignore #2 lines.
            elif line[0:2] == '#2':
                donothing = 1
            #split .bnx files are now single color. Write down the #QX1 lines in each file.
            elif line[0:19] == '# Quality Score QX1':
                gfile.write(line)
                rfile.write(line)
            #split .bnx files are now single color. ignore the #QX2 lines in each file.
            elif line[0:19] == '# Quality Score QX2':
                donothing = 1
            #else... write down all other # lines in the split .bnx files. These lines are not color dependent.    
            else:
                gfile.write(line)
                rfile.write(line)                           

        #Split molecule lines (0,1,2,QX11,QX21 etc.) per molecule based on color.
        else:

            if line[0] == '0':
                #store line in line0 if first entry is 0.
                line0 = line

            if line[0] == '1':
                #if line is 1, determine the number of labels. -2 because ignore first value of 1 and last value molecule length.
                word = line.split('\t')
                numLabels = len(word) - 2                   

                #in line0 there is a number of labels value green+red. Want to find number of green or red labels only when writing new 0 line for split .bnx
                word = (line0.strip()).split('\t')
                word[5] = numLabels #edit total number of labels in line 0 with new value for green or red.
                
                #write new 0 line with edited green or red number of labels based on color of 1 line.
                if colorHash['1'] == 'Green':
                    gfile.write('0')
                    for i in range(1,len(word)):
                        gfile.write('\t'+str(word[i]))
                    gfile.write('\n')
                    
                    #write un-edited 1 line to split file.
                    gfile.write(line)
                elif colorHash['1'] == 'Red':
                    rfile.write('0')
                    for i in range(1,len(word)):
                        rfile.write('\t'+str(word[i]))
                    rfile.write('\n')

                    #write un-edited 1 line to split file.
                    rfile.write(line)                    

            elif line[0] == '2':    #Same logic as 1 line
                #if line is 2, determine the number of labels. -2 because ignore first value of 1 and last value molecule length.
                word = line.split('\t')
                numLabels = len(word) - 2                   

                #in line0 there is a number of labels value green+red. Want to find number of green or red labels only when writing new 0 line for split .bnx
                word = (line0.strip()).split('\t')
                word[5] = numLabels

                if colorHash['2'] == 'Green':
                    gfile.write('0')
                    for i in range(1,len(word)):
                        gfile.write('\t'+str(word[i]))
                    gfile.write('\n')

                    #write un-edited 2 line to split file.
                    gfile.write('1'+line[1:])
                elif colorHash['2'] == 'Red':
                    rfile.write('0')
                    for i in range(1,len(word)):
                        rfile.write('\t'+str(word[i]))
                    rfile.write('\n')    

                    #write un-edited 2 line to split file.                             
                    rfile.write('1'+line[1:])

            #if QX1- code, write to appropriate new file.
            if line[0:4] == 'QX11':
                if colorHash['1'] == 'Green':
                    gfile.write(line)
                elif colorHash['1'] == 'Red':
                    rfile.write(line)

            elif line[0:4] == 'QX12':
                if colorHash['1'] == 'Green':
                    gfile.write(line)
                elif colorHash['1'] == 'Red':
                    rfile.write(line)
 
            #if QX2- code, rewrite as a QX1- code to appropriate file.
            elif line[0:4] == 'QX21':
                if colorHash['2'] == 'Green':
                    gfile.write('QX11'+line[4:])
                elif colorHash['2'] == 'Red':
                    rfile.write('QX11'+line[4:])

            elif line[0:4] == 'QX22':
                if colorHash['2'] == 'Green':
                    gfile.write('QX12'+line[4:])
                elif colorHash['2'] == 'Red':
                    rfile.write('QX12'+line[4:])

    infile.close()
    gfile.close()
    rfile.close()


parser = argparse.ArgumentParser(description='Splits a 2 color .bnx file into a separate green and red .bnx file.')
parser.add_argument("-b", "--bnxFile", help="2 color bnx file name, e.g. RawMolecules.bnx",type=str,default='')
parser.add_argument("-p", "--prefix", help="prefix for output green or red bnx file.",type=str,default='some_great_data')
args = parser.parse_args()

#function reads 2 color .bnx file and outputs a green and red .bnx file.
splitBNXFile(args.bnxFile,args.prefix)

