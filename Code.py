import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import csv
import sp

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
    def __init__(self,FileName=None,StarName="",SpectralType=None,SpectralNumber=None,
    LumType=None,WaveLength=[],Flux=[]):
    
        self.StarName = StarName
        self.SpectralType = SpectralType
        self.SpectralNumber = SpectralNumber
        self.LumType = LumType
        
        if(len(WaveLength) != 0 and len(Flux) != 0):
            self.WaveLength = WaveLength #armstrong
            self.Flux = Flux #erg/cm^2/s/armstrong
        else:
            if(FileName != None):
                self.WaveLength, self.Flux = ReadFile(FileName) 
                #Gets wavelength,flux from file
        
        
        self.GetColors() #Gets the colors of the stars
        
        
    def PlotFlux(self,Xmin=0,Xmax=None,Xlabel=r"Wavelength($\AA$)",Ylabel=r"Flux($erg/cm^2/s/Hz$)"): 

        plt.style.use("seaborn-darkgrid")
        plt.figure(figsize=[8,4.5])
        
        if(Xmax==None): #If no Xmax is given, Use maximum of wavelength
            Xmax = np.amax(WaveLength)

        plt.plot(self.WaveLength,self.Flux)

        if(self.LumType == "I"):
            name = "giant"
        else:
            name = "main sequence star"
            
        plt.title("Flux of a {}{} {} as a function of Wavelength".format(
        self.SpectralType,self.SpectralNumber,name))
        
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

def ReadFile(FileName,Dtype='float'): #Reads a file containing Flux(lambda)
    TotalName = prefix+FileName

    WaveLength, Flux = np.loadtxt(TotalName, delimiter=",",unpack=True,dtype=Dtype)
    return WaveLength, Flux


def EnlistStars(FileNames,StarTypes): #Creates a list if star objects for all stars in files
    
    StarList = np.empty(len(FileNames),dtype=Star)
    for i in range(len(FileNames)):
        LumType = StarTypes[i][3] #Luminosity class
        SpectralNumber = StarTypes[i][2] #ex. the 5 in A5
        SpectralType = StarTypes[i][1] #ex. the A in A5
        StarName = FileNames[i].strip(".dat.fix") #remove filetype
        StarList[i] = Star(FileNames[i],StarName,SpectralType,SpectralNumber,LumType)
        
    return StarList


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

def Integrate(x,y,f,Xmin=None,Xmax=None):  #calculates he integral of f(x)dx from Xmin to Xmax with trapezoids

    if(not(Xmin == None and Xmax == None)):
        Interval = (x>Xmin)*(x<Xmax)
        x = x[Interval]
        y = y[Interval]

    Result = 0
    for i in range(len(x)-1):
        Result += 0.5*(x[i+1]-x[i])*(y[i]+y[i+1])*f((x[i]+x[i+1])/2) #surface ofa trapezoid times f(x)

    return Result


#Wavelength in armstrong, Flux in erg/s/cm^/Hz
def ConvertToWaveLength(WaveLength,Flux): #Converts unit of Flux erg/s/cm^2/HZ -> erg/s/cm^2/armstrong
    
    DeltaLambda = (WaveLength/1E10)**2 / c #m
    OneArmstrong = 1E-10 #m
    return Flux * (OneArmstrong / DeltaLambda) #Now Flux is in erg/s/cm^2/armstrong


def f1(WaveLength): #Required to put into the inegratge function
    return WaveLength;
    

def CalculateC(WaveVega,FluxVega,LambdaFilter,DeltaFilter): #Calculates the c constant for a filter
    
    FilterMin = LambdaFilter - DeltaFilter
    FilterMax = LambdaFilter + DeltaFilter
    Numerator = Integrate(WaveVega,FluxVega,f1,FilterMin,FilterMax)
    Denominator = 2*LambdaFilter*DeltaFilter #We calculated this analytically
    C = 2.5*np.log10(Numerator/Denominator)
    return C

def Blackbody(WaveLengthOld,T): #
    WaveLength = WaveLengthOld.copy()
    WaveLength *= 1E-10 #From armstrong to m

    return 1E-7 * (((2*h*c**2)/(WaveLength**5))/np.expm1((h*c)/(WaveLength*Kb*T)))
    #Returns flux in erg/cm^2/s/armstrong


def CalcBBColor(Tlist): #Calculates U-B and B-V for blackbody with temperature T
    WL = np.linspace(2999,6768,1000) #armstrong
    U_B = []
    B_V = []
    
    for T in Tlist:
        flux = Blackbody(WL,T) #erg/cm^2/s/armstrong
        BBstar = Star(WaveLength = WL,Flux=flux)
        U_B.append(BBstar.U_B)
        B_V.append(BBstar.B_V)
    
    return U_B,B_V

def T_From_BV(Tlist,BVlist,BV0list): #Calculates Temperature of blackbody with certain B-V color
    TOutputList = []
    
    for BV0 in BV0list:
        i = np.argmin(np.abs(BVlist-BV0)) #Calculates index where BV is closest to B0
        TOutputList.append(Tlist[i])
    return TOutputList

def PlotDividedFlux():

    plt.style.use('seaborn')


    Tlist = np.linspace(2000,30000,500)
    bbU_B, bbB_V = CalcBBColor(Tlist) #List of colors for blackbody's with different T


    BV0List = [] #B-V colours for A0, B0, F0 and G0 stars
    IndexList = [0,18,33,50] #List of indexes of A0,B0,F0,G0 stars
    for i in IndexList:
        BV0List.append(StarList[i].B_V)
    
    #Temperatures of blackbody's with B-V's in BV0List
    T0List = T_From_BV(Tlist,bbB_V,BV0List)
    
    for i in range(len(T0List)): #Plots star flux divided by bb flux 
            plt.figure(figsize=(8,5))
            
            
            ColorList = ["green","blue","orange"]
            FilterNameList = ["U","B","V"]
            for k in range(len(LambdaList)): #Plots vertical Filter lines
                plt.axvline(x=LambdaList[k]-DeltaList[k],linestyle="--",color=ColorList[k],label=FilterNameList[k])
                plt.axvline(x=LambdaList[k]+DeltaList[k],linestyle="--",color=ColorList[k])

            
            star = StarList[IndexList[i]]
            T0 = T0List[i]
            
            bbFlux = Blackbody(star.WaveLength,T0) #Calculates flux of bb
            DivisedFlux = star.Flux / bbFlux
            
            plt.plot(star.WaveLength,DivisedFlux)
            #plt.plot(star.WaveLength,bbFlux/np.amax(bbFlux))
            #plt.plot(star.WaveLength,star.Flux/np.amax(star.Flux))
            
            props = dict(boxstyle='round', facecolor='lightgrey', alpha=1) #
            text = "Tbb = " +str(np.round(T0,0)) #String with used constants
            plt.text(0.7*np.amax(star.WaveLength),0.95*np.amax(DivisedFlux),text,fontsize=11,bbox=props)
                  
            plt.title("Flux of " +str(star.SpectralType)+ "0 star divided by flux of blackbody with same B-V")
            plt.xlabel(r"WaveLength($\AA$)")
            plt.ylabel(r"$F_{star}$ / $F_{bb}$") 
            
            plt.xlim(xmin=2000,xmax=7000)

            plt.legend()

def SelectTemp(T,Tmin,Tmax): #Selects all temperatures with a range
    Tnew = T.copy()
    return Tnew[(Tnew>Tmin) * (Tnew<Tmax)]
    
def DoubleSort(List1,List2):
    
    Copy1 = List1.copy()
    Copy2 = List2.copy()
    
    NewList1 = []
    NewList2 = []
    
    for i  in range(len(List1)-1):
        index = np.argmin(Copy1)
        NewList1.append(Copy1[index])
        NewList2.append(Copy2[index])
        del Copy1[index]
        del Copy2[index]
        
    return NewList1, NewList2


def PlotHertzsprung(LumType): #Plots Hertzsprung russel diagram of stars in StarList

    plt.style.use('dark_background')
    plt.figure()
    

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


    #Plot blackbody colors
    Tlist = np.linspace(2000,30000,2000)
    
    for Spec in SpecDict:
        SelT = SelectTemp(Tlist,SpecDict[Spec].Tmin,SpecDict[Spec].Tmax)
        bbU_B, bbB_V = CalcBBColor(SelT)
        plt.plot(bbU_B, bbB_V,color=SpecDict[Spec].Color)
        
    plt.plot(0,-100,label="Black body")    
    
    #--Indicator--#
    U_BTotal,B_VTotal = DoubleSort(U_BTotal,B_VTotal)
    FitData = np.polyfit(U_BTotal,B_VTotal,6)
    PolyFit = np.poly1d(FitData)
    Xnew = np.linspace(-0.8,1.6,2000)
    
    plt.plot(Xnew,PolyFit(Xnew),linewidth=18,alpha=0.2,color="white")
    #--Indicator--#
    
    #Needed for Title 
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
        
        
def GetStarData(FileName,ExTime,M0,Gain,Radius): #Gets the data from stars in an image
        
    #ExTime = Exposure time
    #M0 = Magnitude zeropoint
    #Gain = photons per ADU
    #Radius = assumed radius of stars in pixels
    
    hdul = fits.open(prefix+FileName) 
    Image = hdul[0].data #Flux for each pixel

    FoundStars = sp.find(Image,1000,5) #Data from the stars in the V image       
    XList = FoundStars[0] #List of x positions for all stars
    YList = FoundStars[1] #List of y positions for all stars
    
    #Gets the flux in ADU for all stars    
    ADU,ADUerr,Sky,SkyErr = sp.aper(Image,XList,YList,Gain,[Radius],[15,25],exact=True,flux=True)
    Flux = ADU / ExTime #ADU / s
    M = M0 -2.5 * np.log10(Flux) #Magnitude
        
    return Image,XList,YList,ADU,ADUerr,M,Sky,SkyErr


def SNRatio(Nsource,Nsky,g,RON,R): #Calculates signal to noise ratio
    k = np.pi * R**2 #Amount of pixels
    Numerator = g*Nsource
    Denominator = np.sqrt(g*Nsource+k*(g*Nsky+RON**2))
    return Numerator / Denominator
    
def Merror(SN):
    return (2.5/np.log(10))/SN

#----------Functions----------#


#----------GlobalVariables----------#

prefix = "C:/Users/Jasper/Desktop/Uni/ALOP/Data/"

#---constants---#
c = 299792458 #m/s
h = 6.626E-34 #Js
Kb = 1.38064852E-23 #J/k
#---constants---#

#---Filters---#
#Ranges of different filters
LambdaU = 3659 #armstrong
DeltaU = 660

LambdaB = 4582
DeltaB = 940

LambdaV = 5448
DeltaV = 880

LambdaList = [LambdaU,LambdaB,LambdaV]
DeltaList = [DeltaU,DeltaB,DeltaV]
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

#print(CU,CB,CV)
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

#---Images---#

R = 5
RON = 250
Gain = 5

Vimage,Xstar,Ystar,V_ADU,V_ADUerr,MV,VSky,VSkyErr = GetStarData("exposure_V.fits",10,17.40,Gain,R)
Bimage,Xstar,Ystar,B_ADU,B_ADUerr,MB,BSky,BSkyErr = GetStarData("exposure_B.fits",30,16.43,Gain,R)

VSNratio = SNRatio(np.transpose(V_ADU),VSky,Gain,RON,R)
BSNratio = SNRatio(np.transpose(B_ADU),BSky,Gain,RON,R)
MV_err = np.round(Merror(VSNratio),2)
MB_err = np.round(Merror(BSNratio),2)

DataTable = np.zeros((5,7))
DataTable[:,0] = np.round(Xstar,0) #X
DataTable[:,1] = np.round(Ystar,0) #Y
DataTable[:,2] = np.round(np.transpose(V_ADU/10)) #V Flux
DataTable[:,3] = np.round(np.transpose(B_ADU/30)) #B Flux
DataTable[:,4] = np.round(VSNratio) #V SNratio
DataTable[:,5] = np.round(BSNratio) #B SNratio
DataTable[:,6] = np.round(np.transpose(MV),3) #V magnitude
#print(np.round(np.transpose(MB),3)) #V magnitude
#print(np.round(np.transpose(MB)-np.transpose(MV),5)) #V magnitude

print("\n")
#print(DataTable)
print(MV_err)
print(MB_err)
#print(Ystar)
#print(np.transpose(VSNratio))

#print("\n SNRatio: "+str(VSNratio))

#---Images---#

#----------GlobalVariables----------#


#----------Main----------#

if(False): #Plot Flux(Wavelength) of measured star
    WaveLength, Flux = ReadFile("Bb11.dat.fix") #,3000,9000,"Flux of STARNAME as a function of wavelength"
    PlotFlux(WaveLength,Flux,3000,9000,"Flux of STARNAME as a function of wavelength")

if(False):
    StarList[3].PlotFlux(Xmin=3000,Xmax=9000)
    StarList[50].PlotFlux(Xmin=3000,Xmax=9000)
    StarList[35].PlotFlux(Xmin=3000,Xmax=9000)


if(False): #Plot Flux(Wavelength) of Vega
    PlotFlux(WaveVega,FluxVega,0,12000,"Flux of Vega as Function of wavelength")
    
if(False): 
    Bb11 = Star("Bb11.dat.fix","Bb11")

if(True):
    PlotHertzsprung("V")
    #PlotHertzsprung("I")
if(False):
    PlotDividedFlux()

if(False): #Plots V-band image of the CCD exposure
    plt.figure()
    
    plt.title("Image of 5 stars in the visible band")
    
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.imshow(np.log(Vimage), cmap='gray', interpolation='nearest')
#----------Main----------#

#----------End----------#

#for i in range(len(StarList)):
#    print(str(i)+": ")
#    print(str(StarList[i].SpectralType)+str(StarList[i].SpectralNumber)+"\n")






plt.show()

#----------End----------#
