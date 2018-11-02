import numpy as np
import matplotlib.pyplot as plt
import csv


#----------SpectralClassObject----------#

class SpectralClass:
    
    def __init__(self,Symbol,Tmin,Tmax,R,G,B):
        self.Symbol = Symbol
        self.Color = self.ColorRGB(R,G,B)
        self.Tmin = Tmin
        self.Tmax = Tmax

    def ColorRGB(self,R,G,B):
        return (R/255,G/255,B/255)

#----------SpectralTypeObject----------#


#----------StarObject----------#

class Star:
    def __init__(self,FileName,StarName="",SpectralType=None,SpectralNumber=None,LumType=None):
        self.StarName = StarName
        self.SpectralType = SpectralType
        self.LumType = LumType
        self.WaveLength, self.Flux = ReadFile(FileName) #Gets wavelength,flux from file
        self.GetColors() #Gets the colors of the stars
        
        
    def PlotFlux(self,Xmin=0,Xmax=None,Xlabel=r"Wavelength($\AA$)",Ylabel=r"Flux($erg/cm^2/s/Hz$)"): 

        plt.figure()    
    
        if(Xmax==None): #If no Xmax is given, Use maximum of wavelength
            Xmax = np.amax(WaveLength)

        plt.plot(self.WaveLength,self.Flux)

        plt.title("Flux of star"+self.StarName+"as a function of Wavelength")
        plt.xlabel(r"Wavelength($\AA$)")
        plt.ylabel(r"Flux($erg/cm^2/s/\AA$)")

        plt.ylim(ymin=0)
        plt.xlim(xmin=Xmin,xmax=Xmax)
    
    def CalculateColor(self,LambdaFilter,DeltaFilter,C):
    
        FilterMin = LambdaFilter - DeltaFilter
        FilterMax = LambdaFilter + DeltaFilter   
        Numerator = Integrate(self.WaveLength,self.Flux,f1,FilterMin,FilterMax)
        Denominator = 2*LambdaFilter*DeltaFilter
        
        Color = -2.5*np.log10(Numerator/Denominator) + C
        return Color

    def GetColors(self):
        self.MU = self.CalculateColor(LambdaU,DeltaU,CU)
        self.MB = self.CalculateColor(LambdaB,DeltaB,CB)
        self.MV = self.CalculateColor(LambdaV,DeltaV,CV)
        self.U_B = self.MU - self.MB
        self.B_V = self.MB - self.MV

#----------StarObject----------#

#----------Functions----------#
def Blackbody(WaveLength,T): #
    WaveLength *= 1E-10 #From armstrong to m
    return ((2*h*c**2)/(WaveLength**5))/np.expm1((h*c)/(WaveLength*Kb*T))
    #Returns flux in w/m^2

def PlotHertzsprung(LumType):

    plt.figure()
    plt.style.use('dark_background')

    U_BTotal = []
    B_VTotal = []

    for Spec in SpecDict:
        U_B = []
        B_V = []
        
        for Star in StarList:
            if(Star.SpectralType == Spec and Star.LumType == LumType):
                U_B.append(Star.U_B)
                B_V.append(Star.B_V)
                U_BTotal.append(Star.U_B)
                B_VTotal.append(Star.B_V)
                
        Color = SpecDict[Spec].Color
        plt.scatter(U_B,B_V,color=Color,s=5,label=Spec+" Type")

        
    if(LumType == "I"):
        startype = "Giants"
    else:
        startype = "Main sequence stars"
    
    plt.title("Hertzsprung-Russel diagram for several "+startype)
    plt.xlabel("U-B")
    plt.ylabel("B-V")
        
    plt.ylim(ymin=1.1*np.amax(B_VTotal),ymax=1.1*np.amin(B_VTotal))
    
    plt.legend()    
    plt.show()
        


def EnlistStars(FileNames,StarTypes): #Creates a list if star objects for all stars in files
    
    StarList = np.empty(len(FileNames),dtype=Star)
    for i in range(len(FileNames)):
        LumType = StarTypes[i][3] #Luminosity class
        SpectralNumber = StarTypes[i][2] #ex. the 5 in A5
        SpectralType = StarTypes[i][1] #ex. the A in A5
        StarName = FileNames[i].strip(".dat.fix") #remove filetype
        StarList[i] = Star(FileNames[i],StarName,SpectralType,SpectralNumber,LumType)
        
    return StarList
    


def FindLumClass(TypeList): #Finds the type from the typeList of steller classification
    LumList = np.empty(len(TypeList),dtype="str")
    for i in range(len(TypeList)):
        LumList[i] = TypeList[i][3]
    return LumList

#Wavelength in armstrong, Flux in erg/s/cm^/Hz
def ConvertToWaveLength(WaveLength,Flux): #Converts unit of Flux erg/s/cm^2/HZ -> erg/s/cm^2/armstrong
    
    DeltaLambda = (WaveLength/1E10)**2 / c #m
    OneArmstrong = 1E-10 #m
    return Flux * (OneArmstrong / DeltaLambda) #Now Flux is in erg/s/cm^2/armstrong


def Integrate(x,y,f,Xmin=None,Xmax=None):  #calculates he integral of f(x)dx from Xmin to Xmax with trapezoids

    if(not(Xmin == None and Xmax == None)):
        Interval = (x>Xmin)*(x<Xmax)
        x = x[Interval]
        y = y[Interval]

    Result = 0
    for i in range(len(x)-1):
        Result += 0.5*(x[i+1]-x[i])*(y[i]+y[i+1])*f((x[i]+x[i+1])/2) #surface ofa trapezoid times f(x)

    return Result


def ReadFile(FileName,Dtype='float'): #Reads a file containing Flux(lambda)
    TotalName = prefix+FileName

    WaveLength, Flux = np.loadtxt(TotalName, delimiter=",",unpack=True,dtype=Dtype)
    return WaveLength, Flux
    

#Plots the flux(wavelength) of a star (NOT DONE) 
def PlotFlux(WaveLength,Flux,Xmin=0,Xmax=None,Title="",Xlabel=r"Wavelength($\AA$)",Ylabel=r"Flux($erg/cm^2/s/Hz$)"): 

    plt.figure()    
    
    if(Xmax==None): #If no Xmax is given, Use maximum of wavelength
        Xmax = np.amax(WaveLength)

    plt.plot(WaveLength,Flux)

    plt.title(Title)
    plt.xlabel(r"Wavelength($\AA$)")
    plt.ylabel(r"Flux($erg/cm^2/s/\AA$)")


    plt.xlim(xmin=Xmin,xmax=Xmax)
    plt.ylim(ymin=0)


def f1(WaveLength): #Required to put into the inegratge function
    return WaveLength;
    

def CalculateC(WaveVega,FluxVega,LambdaFilter,DeltaFilter): #Calculates the c constant for a filter
    
    FilterMin = LambdaFilter - DeltaFilter
    FilterMax = LambdaFilter + DeltaFilter
    Numerator = Integrate(WaveVega,FluxVega,f1,FilterMin,FilterMax)
    Denominator = 2*LambdaFilter*DeltaFilter #We calculated this analytically
    C = 2.5*np.log10(Numerator/Denominator)
    return C


def Blackbody(WaveLength,T): #
    WaveLength *= 1E-10 #From armstrong to m
    return ((2*h*c**2)/(WaveLength**5))/np.expm1((h*c)/(WaveLength*Kb*T))
    #Returns flux in w/m^2
    
#----------Functions----------#


#----------GlobalVariables----------#

prefix = "C:/Users/Jasper/Desktop/Uni/ALOP/Data/"
c = 299792458 #m/s
h = 6.626E-34 #Js

#---Filters---#
LambdaU = 3659 #armstrong
DeltaU = 660

LambdaB = 4582
DeltaB = 940

LambdaV = 5448
DeltaV = 880
#---Filters---#

#---Vega---#
WaveVega, FluxVegaHZ = ReadFile("vega.fnu.csv")#,0,5000
WaveVega *= 10 #Convert from nm to armstrong
FluxVega = ConvertToWaveLength(WaveVega,FluxVegaHZ)
#---Vega---#

#---Cconstants---#
CU = CalculateC(WaveVega,FluxVega,LambdaU,DeltaU)
CB = CalculateC(WaveVega,FluxVega,LambdaB,DeltaB)
CV = CalculateC(WaveVega,FluxVega,LambdaV,DeltaV)

print(CU,CV,CB)

#---Cconstants---#

#---StarList---#
FileNames,StarTypes = ReadFile("spectra.list.csv","str")
StarList = EnlistStars(FileNames,StarTypes)
#---StarList---#

#---SpectralClasses---#


SpecO = SpectralClass("O",30000,1E10,125,151,255)
SpecB = SpectralClass("B",10000,30000,170,191,255)
SpecA = SpectralClass("A",7500,10000,213,244,255)
SpecF = SpectralClass("F",6000,7500,249,245,255)
SpecG = SpectralClass("G",5200,6000,255,237,227)
SpecK = SpectralClass("K",3700,5200,255,218,181)
SpecM = SpectralClass("M",2400,3700,255,181,108)

SpecDict = {
    "O": SpecO,
    "B": SpecB,
    "A": SpecA,
    "F": SpecF,
    "G": SpecG,
    "K": SpecK,
    "M": SpecM,
    
}

#---SpectralClasses---#

#----------GlobalVariables----------#





#----------Main----------#

if(False): #Plot Flux(Wavelength) of measured star
    WaveLength, Flux = ReadFile("Bb11.dat.fix") #,3000,9000,"Flux of STARNAME as a function of wavelength"
    PlotFlux(WaveLength,Flux,3000,9000,"Flux of STARNAME as a function of wavelength")


if(True): #Plot Flux(Wavelength) of Vega
    PlotFlux(WaveVega,FluxVega,0,12000,"Flux of Vega as Function of wavelength",)
    
if(True): 
    Bb11 = Star("Bb11.dat.fix","Bb11")

if(True):
    PlotHertzsprung("V")
    PlotHertzsprung("I")

#----------Main----------#

#----------End----------#

plt.show()

#----------End----------#