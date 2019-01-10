import random
import numpy as np
import matplotlib.pyplot as plt
import itertools

def CollisionFinder(MCLine,MCDir,LinePosition,MCBoundary,BoundaryLimits):
    '''
    COLLISION FINDER
    Function which calculates the collision POSITION given: the line gradient (M) and constant (C), the line\'s previous position
    The line\'s direction (To differentiate between the two direction-based solutions), then the boundary gradient and constant
    and the limits of the boundary.

    MCLine: a tuple of two integer values (gradient,constant), defining the line of motion of the particle
    MCDir: the x-direction of the line, takes one of two values: +1 for increasing x, -1 for decreasing x
    LinePosition: (x,y) tuple of the previous position of the particle (i.e usually the last wall position of the collision)
    MCboundary: a tuple of two integer values (gradient,constant), defining the line of the boundary
    BoundaryLimits: two lists nested within a list, each containing a set of [xmin,xmax],[ymin,ymax] values for the boundary conditions

    Return Values:
    Ideal output: an (x,y) tuple of the new collision position
    Else: returns a string stating why no collision was made

    To use this function in a closed system iterate over all boundaries
    '''
    if MCLine[0]==MCBoundary[0]:
        #print('No Collision - Same Grad')
        return 'No Collision - Same Grad'

    XSolution = (MCBoundary[1]-MCLine[1])/(MCLine[0]-MCBoundary[0])
    YSolution = XSolution*MCLine[0]+MCLine[1]
    if BoundaryLimits==None:
        return [XSolution,YSolution,'Special']
    if (LinePosition[0]<min(BoundaryLimits[0]) and MCDir==-1) or (LinePosition[0]>max(BoundaryLimits[0]) and MCDir==1):
        #print('No Collision - Wrong Direction')
        return 'No Collision - Wrong Direction'
    if XSolution>max(BoundaryLimits[0]) or XSolution<min(BoundaryLimits[0]) or YSolution>max(BoundaryLimits[1]) or YSolution<min(BoundaryLimits[1]):
        #print('No Collision - Out of Boundary Limits')
        return 'No Collision - Out of Boundary Limits'
    if (XSolution>LinePosition[0] and MCDir==-1) or (XSolution<LinePosition[0] and MCDir==1):
        #print('No Collision - Wrong Direction 2')
        return 'Non Collision - Wrong Direction 2'
    return [XSolution,YSolution]

def Reflector(ParticleAngle,WallAngle,Randomness=0):

    '''
    Function that returns the reflected angle of an incoming line and a wall.
    We define 0 degrees as the positive x direction, increasing anticlockwise. All angles are between 0 and 359.999...

    ParticleAngle: The incoming angle of the particle in degrees
    WallAngle: The angle of the wall, since this does not have a direction the wall angle can be either direction along its line
               e.g the 45 degree wall, given by y=x line, will return the same reflection as a 225 degree wall
    Randomness: Optional arg which randomises the return angle. This will be limited by the wall angle so that particle doesn't ever go through
                the wall.

    Return Values:
    Only output: single integer of the new reflected outgoing angle.
    '''

    WallAngle = (WallAngle+360)%180
    OutAngle = ((2*WallAngle - ParticleAngle)+360)%360
    LimitAngle = min(abs(WallAngle-OutAngle),abs(180+WallAngle-OutAngle),abs(360+WallAngle-OutAngle))
    Randomness = min(abs(Randomness),round(LimitAngle/1.5))
    #print('ParticleAngle: '+str(ParticleAngle)+ '   WallAngle: '+str(WallAngle)+ '    OutAngle:'+str(OutAngle))
    ReturnedAngle = (OutAngle + random.randint(-Randomness,Randomness))%360
    while ReturnedAngle%90==0:
        ReturnedAngle = (OutAngle + random.randint(-Randomness,Randomness))%360
    return ReturnedAngle

def ConvertDegToMC(Degrees):

    '''
    Function that takes a degree value and tranforms it into line gradient information and a XDir to differentiate between two possible directions

    Degrees: integer input between 0 and 360

    Return values:
    Only return: A combination of (gradient,xdirection)
    Note: No constant returned for line information
    '''

    XDir = -1 if (Degrees <= 270 and Degrees >90) else  1
    return np.tan(Degrees*np.pi/180),XDir

def ConvertMCToDeg(MLine,XDir):

    '''
    Function that takes a line gradient and direction and transforms it into an angle in degrees

    MLine: Gradient information about the line
    XDir: XDirection, takes +1 if increasing x, else -1

    Return values:
    Only return: A integer value of degrees
    '''
    Additional = 180 if (XDir == -1) else 0
    Additional = 360 if (not Additional and np.sign(MLine)==-1) else Additional
    return ((np.arctan(MLine)*180/np.pi) + Additional)%360

def DistanceChooser(StartPosition,HitPositions,HitNumbers):

    '''
    Given all of the collision points, this function will return the closest one to the origin, this helps create a more
    realistic situation to avoid a particle travelling through a wall and colliding with the one behind it instead

    StartPosition: a tuple of (x,y) values for the previous position of the particle
    HitPositions: a list of all (x,y) collisions found
    HitNumbers: a list of integers of the collisional walls\'s indeces within the total boundary limits, this helps determine which
                wall should be ignored over the next iteration of collisions after colliding with a wall

    Return values:
    Ideal output: the position (x,y) of the nearest collision and the index of the wall in the total boundaries list
    Else: if the HitPositions is an empty list then no collisions have been found and an error is raised
    '''

    if len(HitPositions)==0:
        print('SOMETHING WRONG NOTHING HIT')
        k=input()

    if len(HitPositions)==1:
        if len(HitPositions[0])==3:
            return HitPositions[0][:2],'Break'
        return HitPositions[0],HitNumbers[0]
    func = lambda x:((x[0]-StartPosition[0])**2 + (x[1]-StartPosition[1])**2)**0.5
    MinDistance= min(HitPositions,key=func)
    #print('HN'+str(HitNumbers))
    #print(HitNumbers[HitPositions.index(MinDistance)])
    return MinDistance,HitNumbers[HitPositions.index(MinDistance)]

def MToMC(M, Position):
    '''
    Converts a gradient and a position into (gradient,constant) line information

    M: the gradient of the line
    Position: the current (x,y) position of the particle

    Return values:
    Only output: (gradient,constant) integers
    '''
    Constant = Position[1]-M*Position[0]
    return [M,Constant]

def MakeCapillaries(Number,X1):
    #Makes the capillary arrays for how ever many capills we want, woooooo
    FunctCapillaries=[]
    FunctCapillaryBounds=[]
    Step = 6/Number
    for i in range(1,Number):
        Boundary = [0,i*Step-3]
        BoundaryLim = [[X1,X1+20],[i*Step-2.999999,i*Step-3.0000001]]
        FunctCapillaries.append(Boundary)
        FunctCapillaryBounds.append(BoundaryLim)
    print('{} boundaries made! At y={}'.format(Number-1,list(i*Step-3 for i in range(1,Number))))
    return FunctCapillaries,FunctCapillaryBounds

def main(IterationNumber,TotalBoundaries,TotalLims):
    '''
    Main Iterator, does shit, not sure how
    '''
    #CONDITIONS
    RandomAngle = 20
    Start = [150,0]
    StartAng = random.randint(0,359)
    StartM,StartDir = ConvertDegToMC(StartAng)
    StartMC = MToMC(StartM,Start)

    #print(TotalBoundaries)

    ###   Initiating data stores   ####
    XYPlotHold=[[],[]]
    XYPlotHold[0].append(Start[0])
    XYPlotHold[1].append(Start[1])
    CurrentPosition=Start
    CurrentMC = StartMC
    CurrentDir = StartDir
    CurrentAng = ConvertMCToDeg(CurrentMC[0],CurrentDir)
    IgnoreWall=-1

    ###   Iterations    ###
    for _ in range(IterationNumber):
        PossibleHits=[]
        HitNumber=[]
        for NIndex,(Bound,BLims) in enumerate(zip(TotalBoundaries,TotalLims)):
            #print('Iteration: '+ str(_) + '    NIndex: '+str(NIndex))
            if NIndex==IgnoreWall:
                continue
            #Iterates over all possible hits, must choose nearest
            Coll = CollisionFinder(CurrentMC,CurrentDir,CurrentPosition,Bound,BLims)
            if type(Coll)!=str:
                if len(Coll)==2:
                    PossibleHits.append(Coll)
                    HitNumber.append(NIndex)
                if len(Coll)==3:
                    PossibleHits.append(Coll)
                    HitNumber.append(NIndex)
        CurrentPosition,IgnoreWall=DistanceChooser(CurrentPosition,PossibleHits,HitNumber)
        if IgnoreWall=='Break':
            XYPlotHold[0].append(CurrentPosition[0])
            XYPlotHold[1].append(CurrentPosition[1])
            break
        CurrentAng=Reflector(CurrentAng,ConvertMCToDeg(TotalBoundaries[IgnoreWall][0],1),RandomAngle) #Randomness for set above
        CurrentM,CurrentDir=ConvertDegToMC(CurrentAng)
        CurrentMC=MToMC(CurrentM,CurrentPosition)
        XYPlotHold[0].append(CurrentPosition[0])
        XYPlotHold[1].append(CurrentPosition[1])

        #print('POS:'+str(CurrentPosition)+'   ANG:'+str(CurrentAng)+'    MC:' + str(CurrentMC) + '    DIR: '+str(CurrentDir))

    # if abs(XYPlotHold[1][-1])<2000 and _!=(IterationNumber-1):
    #     plt.plot(XYPlotHold[0],XYPlotHold[1])
    return XYPlotHold

### CONDITIONS ###
X1,X2=120,160
Y1,Y2=-3,3

Boundary1 = [1e10,-X2*1e10]
BoundaryLims1 = [[X2-0.000001,X2+0.000001],[Y1-0.000001,Y2+0.0000001]]

Boundary2 = [0,Y2]
BoundaryLims2 = [[X1,X2+0.00001],[Y1-0.0000001,Y1+0.000001]]

Boundary3 = [0,Y1]
BoundaryLims3 = [[X1,X2+0.00001],[Y2-0.0000001,Y2+0.000001]]

###For A Funnel on the end###
# BoundaryF1 = [1e10,-60e10]
# BoundaryLimsF1 = [[59.99999,60.00001],[1,3.0001]]
#
# BoundaryF2 = [1e10,-60e10]
# BoundaryLimsF2 = [[59.99999,60.00001],[-1,-3.0001]]

### For Capillaries
BoundaryC1=[0,2]
BoundaryLimsC1=[[60,80],[1.999999,2.0000001]]
BoundaryC2=[0,1]
BoundaryLimsC2=[[60,80],[0.99999,1.0000001]]
BoundaryC3=[0,0]
BoundaryLimsC3=[[60,80],[-0.000001,0.0000001]]
BoundaryC4=[0,-1]
BoundaryLimsC4=[[60,80],[-0.999999,-1.00000001]]
BoundaryC5=[0,-2]
BoundaryLimsC5=[[60,80],[-1.999999999,-2.0000001]]
BoundaryCapillaries=[BoundaryC1,BoundaryC2,BoundaryC3,BoundaryC4,BoundaryC5]
BoundaryLimsCapillaries=[BoundaryLimsC1,BoundaryLimsC2,BoundaryLimsC3,BoundaryLimsC4,BoundaryLimsC5]

### For Finer Capillaries (2x)
BoundaryFC1=[0,2.5]
BoundaryLimsFC1=[[60,80],[2.49999999,2.500000001]]
BoundaryFC2=[0,1.5]
BoundaryLimsFC2=[[60,80],[1.49999999,1.500000001]]
BoundaryFC3=[0,0.5]
BoundaryLimsFC3=[[60,80],[0.49999999,0.500000001]]
BoundaryFC4=[0,-0.5]
BoundaryLimsFC4=[[60,80],[-0.49999999,-0.500000001]]
BoundaryFC5=[0,-1.5]
BoundaryLimsFC5=[[60,80],[-1.49999999,-1.500000001]]
BoundaryFC6=[0,2.5]
BoundaryLimsFC6=[[60,80],[-2.49999999,-2.500000001]]
BoundaryFineCapillaries=[BoundaryFC1,BoundaryFC2,BoundaryFC3,BoundaryFC4,BoundaryFC5,BoundaryFC6]
BoundaryLimsFineCapillaries=[BoundaryLimsFC1,BoundaryLimsFC2,BoundaryLimsFC3,BoundaryLimsFC4,BoundaryLimsFC5,BoundaryLimsFC6]

BoundaryEnd = [1e10,0]
#BoundaryLims5 = [[-20.00001,0.00001],[100.00001,99.99999]]
BoundaryLimsEnd = None

FunctCapBounds,FunctCapLims=MakeCapillaries(5,X1)

TotalBoundaries=[Boundary1,Boundary2,Boundary3,BoundaryEnd] + FunctCapBounds # + BoundaryCapillaries + BoundaryFineCapillaries
TotalLims = [BoundaryLims1,BoundaryLims2,BoundaryLims3,BoundaryLimsEnd] + FunctCapLims # + BoundaryLimsCapillaries + BoundaryLimsFineCapillaries



with open('FarCapillCollY.txt',mode='a') as FileWall:
    FileWall.truncate(0)
    YHistHold=[]
    for __ in range(50000):
        if __%500==0:
            print(int(__/500),end='\r')
        YHold=main(200,TotalBoundaries,TotalLims)[1]
        #print(YHold)
        if abs(YHold[-1])<2000 and len(YHold)!=200:
            YHistHold.append(YHold[-1])

            FileWall.write(str(YHold[-1])+'\n')


    plt.axvline(0,color='k')
    plt.plot((60,100,100,60),(-3,-3,3,3),color='k')
    for i in range(-2,3):
        plt.plot([60,80],[i,i],color='k')
    # with open('CapillCollY.txt',mode='a') as filesaveY:
    #     filesaveY.truncate(0)
    #     for YVal in YHistHold:
    #         filesaveY.write(str(YVal)+'\n')
    #print(YHistHold)
    plt.xlabel('X-Position /mm')
    plt.ylabel('Y-Position /mm')
    plt.show()
    plt.hist(YHistHold,200)
    plt.xlabel('Final Y-position /mm')
    StDev=np.std(YHistHold)
    print(StDev)
    #plt.text(max(YHistHold)-400,300,'Standard Deviation: '+str(round(StDev,1)))

    plt.show()
