#!/usr/bin/env python

import urllib
import string
import os
import sys
import LaunchOnCondor
import glob

def ComputeLimits(InputPattern):
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f1"' , '"Gluino300"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f1"' , '"Gluino400"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f1"' , '"Gluino500"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f1"' , '"Gluino600"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f1"' , '"Gluino700"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f1"' , '"Gluino800"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f1"' , '"Gluino900"', 0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f1"', '"Gluino1000"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100_f1"', '"Gluino1100"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1200_f1"', '"Gluino1200"',0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015])

        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f5"' , '"Gluino300"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f5"' , '"Gluino400"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f5"' , '"Gluino500"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f5"' , '"Gluino600"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f5"' , '"Gluino700"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f5"' , '"Gluino800"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f5"' , '"Gluino900"', 0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f5"', '"Gluino1000"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100_f5"', '"Gluino1100"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1200_f5"', '"Gluino1200"',0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015])

        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino300_f10"' , '"Gluino300"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino400_f10"' , '"Gluino400"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino500_f10"' , '"Gluino500"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino600_f10"' , '"Gluino600"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino700_f10"' , '"Gluino700"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino800_f10"' , '"Gluino800"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino900_f10"' , '"Gluino900"', 1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1000_f10"', '"Gluino1000"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1100_f10"', '"Gluino1100"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"ANALYSE"', InputPattern, '"Gluino1200_f10"', '"Gluino1200"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015])

if len(sys.argv)==1:
	print "Please pass in argument a number between 0 and 2"
        print "  0 - Submit the Core of the (TkOnly+TkTOF) Analysis     --> submitting 2x 3 jobs"
        print "  1 - Run the control plot macro                         --> submitting    0 jobs"
        print "  2 - Run the Optimization macro based on best Exp Limit --> submitting 2xSignalPoints jobs"
        print "  3 - Run the exclusion plot macro                       --> submitting    0 jobs"
	sys.exit()

elif sys.argv[1]=='0':	
     print 'ANALYSIS'
     FarmDirectory = "FARM"
     JobName = "HscpAnalysis"
     LaunchOnCondor.Jobs_RunHere = 1
     LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)

     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"190001_191000"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"190001_191000"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191001_192000"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191001_192000"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"192001_193751"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"192001_193751"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"193752_194076"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"193752_194076"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"194077_194300"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"194077_194300"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"194301_194500"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"194301_194500"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"194501_194800"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"194501_194800"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"194701_194900"', 80, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"194701_194900"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"194901_195016"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"194901_195016"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"195099_195378"', 80, 2.1])
     LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"195099_195378"', 80, 2.1])

     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"190645-191100"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191101-191246"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191247"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191248-191300"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191301-191800"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"191801-192000"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"193001-193336"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"193338-193557"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_DATA"' ,'"193558-194076"', 70, 2.1])


     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino300"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino400"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino500"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino600"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino700"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino800"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino900"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino1000"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino1100"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Gluino1200"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau100"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau126"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau156"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau200"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau247"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau308"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau370"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau432"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"GMStau494"',  70, 2.1])

     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop130"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop200"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop300"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop400"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop500"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop600"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop700"',  70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C", '"ANALYSE_SIGNAL"', '"Stop800"',  70, 2.1])

     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"190645-191100"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191101-191246"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191247"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191248-191300"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191301-191800"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"191801-192000"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"193001-193336"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"193338-193557"', 70, 2.1])
     #LaunchOnCondor.SendCluster_Push(["FWLITE", os.getcwd()+"/Analysis_Step234.C",'"ANALYSE_COSMIC"' ,'"193558-194076"', 70, 2.1])

     LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='1':
        print 'PLOTTING'
	os.system('root Analysis_Step5.C++ -l -b -q')

elif sys.argv[1]=='2':
        print 'OPTIMIZATION'
        FarmDirectory = "FARM"
	JobName = "HscpOptmization"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino300_f10"', '"Gluino300"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino400_f10"', '"Gluino400"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino500_f10"', '"Gluino500"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino600_f10"', '"Gluino600"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino700_f10"', '"Gluino700"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino800_f10"', '"Gluino800"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino900_f10"', '"Gluino900"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino1000_f10"', '"Gluino1000"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino1100_f10"', '"Gluino1100"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"/Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Gluino1200_f10"', '"Gluino1200"',1.0 / 0.3029 , 0.0 / 0.4955 , 0.0 / 0.2015 ])

        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop130"'      , '"Stop130"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop200"'      , '"Stop200"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop300"'      , '"Stop300"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop400"'      , '"Stop400"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop500"'      , '"Stop500"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop600"'      , '"Stop600"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop700"'      , '"Stop700"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"Stop800"'      , '"Stop800"', 0.1705 / 0.1705 , 0.4868 / 0.4868 , 0.3427 / 0.3427   ])

        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau100"', '"GMStau100"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau126"', '"GMStau126"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau156"', '"GMStau156"' ,-1, -1, -1  ])
	LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau200"', '"GMStau200"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau247"', '"GMStau247"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau308"', '"GMStau308"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau370"', '"GMStau370"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau432"', '"GMStau432"' ,-1, -1, -1  ])
        LaunchOnCondor.SendCluster_Push(["ROOT", os.getcwd()+"//Analysis_Step6.C", '"FINDCUT"', '"Results/Eta21/PtMin80/"', '"GMStau494"', '"GMStau494"' ,-1, -1, -1  ])

        LaunchOnCondor.SendCluster_Submit()

elif sys.argv[1]=='3':
        print 'LIMITS'
        FarmDirectory = "FARM"
        JobName = "HscpLimits"
        LaunchOnCondor.Jobs_RunHere = 1
        LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)


        ComputeLimits('"Results/Eta21/PtMin80/"')
        LaunchOnCondor.SendCluster_Submit()


elif sys.argv[1]=='4':
        print 'EXCLUSION'
        os.system('sh Analysis_Step6.sh')
else:
	print 'Unknwon case: use an other argument or no argument to get help'



