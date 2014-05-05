# Script Name: Optimized Pit Removal Version 1.5
# 
# Created By:  Stephen Jackson
# Date:        01/11/2013

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
executablename = '\OptimizedPitRemoval.exe'
executablestr = '"' + executablepath + executablename + '"'
arcpy.AddMessage(executablestr)

# Get Input DEM (ASCII)
#
inLyr = arcpy.GetParameterAsText(0)
desc = arcpy.Describe(inLyr)
InputDEM=str(desc.catalogPath)
arcpy.AddMessage("\nInput Elevation file: "+InputDEM)
InputDEMstr = ' -z "' + InputDEM + '"'

# Get Output DEM (ASCII)
OutputDEM = arcpy.GetParameterAsText(1)
Ext = OutputDEM[-4:]
if not (Ext == '.txt'):
    OutputDEM = OutputDEM + '.txt'

arcpy.AddMessage("\nOutput Elevation file: "+OutputDEM)
OutputDEMstr = ' -fel "' + OutputDEM + '"'

# Get Output Pit (ASCII)
#pitFile = arcpy.GetParameterAsText(2)
#if pitFile:
#    Ext = pitFile[-4:]
#    if not (Ext == '.txt'):
#        pitFile = pitFile + '.txt'
#    
#    OuputPitstr = ' -pit "' + pitFile + '"'
#    arcpy.AddMessage("\nOutput Pit Location file: "+pitFile)
#else:
#    arcpy.AddMessage("\nNo Pit Location file selected. Pit locations will not be stored.")

# Get Mode (String)
Mode = arcpy.GetParameterAsText(2)
arcpy.AddMessage("\nOptimization Mode: "+Mode)
if (Mode == 'Minimize Net Elevation Change'):
    Mode = 'bal'
elif (Mode == 'Cut Only'):
    Mode = 'cut'
else:
    Mode = 'mincost'
Modestr = ' -mode ' + Mode

# Get Step Size (double)
StepSize = arcpy.GetParameterAsText(3)
arcpy.AddMessage("\nStep Size: "+StepSize)
StepSizestr = ' -step ' + StepSize

# Construct the command line.  Put quotes around file names in case there are spaces
#if pitFile:
#    cmd = executablestr + InputDEMstr + OutputDEMstr + OuputPitstr + Modestr + StepSizestr
#else:
cmd = executablestr + InputDEMstr + OutputDEMstr + Modestr + StepSizestr
arcpy.AddMessage(cmd)
os.system(cmd)
pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# print the output lines as they happen
while True:
    line = pipe.stdout.readline()
    if not line:
        break
    arcpy.AddMessage(line)
