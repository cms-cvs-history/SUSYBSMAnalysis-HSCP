
#ifndef HSCP_ANALYSIS_GLOBAL
#define HSCP_ANALYSIS_GLOBAL

const double pi = 3.1415926535;

std::string        MODE            = "COMPILE";

int files=9;

double             PtHistoUpperBound   = 3500;
double             MassHistoUpperBound = 2000;
int		   MassNBins           = 200;
double             IPHistoUpperBound   = 200;
double             IPLimit=300;
//double             MinCosmicPt=10;

float              GlobalMaxV3D  =   10;

float              DTRegion      =   0.9;
float              CSCRegion    =   0.9;
float              GlobalMaxDxy  =   15.;
float              GlobalMaxDz   =   12.;
float              CosmicMinDz   =   70.;
float              CosmicMaxDz   =   120.;
float              CosmicMinV3D  =   70.;
float              GlobalMaxDXY  =   10.00;
float              GlobalTkMaxDXY  =   0.2;
float              GlobalMaxChi2 =   5.0;
double             GlobalMinNDOF =   8;
double             GlobalMinNDOFDT  =  6;
double             GlobalMinNDOFCSC =  6;
double             GlobalMaxTOFErr =   0.07;
double             GlobalMaxPterrSq=   9999;
double             GlobalMinPt   =   80.00;
double             GlobalMaxP   =    99999999.;
double             GlobalMinTOF  =   1.0;
float              GlobalMaxEta  =   2.4;
double             MaxDistTrigger=   0.4;
//double             maxSegSep=0.3;
double             minSegEtaSep=0.04;

const int DzRegions=6;
std::string RegionNames[DzRegions]={"Region0","Region1","Region2","Region3","Region4", "Region5"};
std::string LegendNames[DzRegions]={"dz < 6 cm","6 cm < dz < 30 cm","30 cm < dz < 50 cm","50 cm < dz < 70 cm","70 cm < dz < 120 cm", "dz > 120 cm"};


/*
std::string        dEdxM_Label     = "dedxHarm2";

float              GlobalMaxV3D  =   0.50;
float              GlobalMaxChi2 =   5.0;
int                GlobalMinQual =   2;
unsigned int       GlobalMinNOH  =   11;
unsigned int       GlobalMinNOM  =   6;
double             GlobalMinNDOF =   8;
double             GlobalMinNDOFDT  =  6;
double             GlobalMinNDOFCSC =  6;
double             GlobalMaxTOFErr =   0.07;
double             GlobalMaxPterr=   0.25;
double             GlobalMaxTIsol = 50;
double             GlobalMaxEIsol = 0.30;
double             GlobalMinPt   =   50.00;
double             GlobalMinIs   =   0.0;
double             GlobalMinIm   =   3.0;
float              GlobalMaxEta  =  1.5;
*/
#endif
