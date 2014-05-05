# Script Name: CutFillStatistics 1.0
# 
# Created By:  Stephen Jackson
# Date:        01/16/2013

# Import ArcPy site-package and os modules
#
import arcpy
import os
import sys
import time
import string
import subprocess

#set executable program location
executablepath = os.path.dirname(os.path.abspath(__file__))
arcpy.AddMessage(executablepath)
executablename = '\CutFillStatistics.exe'
executablestr = '"' + executablepath + executablename + '"'
arcpy.AddMessage(executablestr)

# Get Original DEM (ASCII)
#
inLyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inLyr)
OrigDEM=str(desc.catalogPath)
arcpy.AddMessage("\nOriginal Elevation file: "+OrigDEM)
OriginalDEMstr = ' -orig "' + OrigDEM + '"'
arcpy.AddMessage(OriginalDEMstr)

# Get Modified DEM (ASCII)
ModDEM = arcpy.GetParameterAsText(1)
arcpy.AddMessage("\Modified Elevation file: "+ModDEM)
ModDEMstr = ' -mod "' + ModDEM + '"'

# Get Output Statistics File (ASCII)
StatFile = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nOutput Statistics file: "+StatFile)
StatFilestr = ' -stat "' + StatFile + '"'


# Construct the command line.  Put quotes around file names in case there are spaces
cmd = executablestr + OriginalDEMstr + ModDEMstr + StatFilestr
arcpy.AddMessage(cmd)
os.system(cmd)
pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# print the output lines as they happen
while True:
    line = pipe.stdout.readline()
    if not line:
        break
    arcpy.AddMessage(line)
