from pcraster import *
from pcraster.framework import *
import csv
import sys 

numberOfTimesExecuted = sys.argv[1]
initialPreyDensity = 0.0
initialPredDensity = 0.0

def printStatus():
    print()
    print("Running the model " + numberOfTimesExecuted +" times")
    print("Prey distribution of:"+str(initialPreyDensity))
    print("Predator distribution of:"+str(initialPredDensity))
    print()

def printToCSV(avg_prey,avg_pred,avg_inf):
    print('Writing to csv now')
    # write to csv 
    with open('csv/predatorInf.csv','a',newline="\n") as file:
        writer  = csv.writer(file, delimiter=';',
                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow([round(initialPreyDensity,1),round(initialPredDensity,1),avg_prey,avg_pred,avg_inf])
    print('Done!')
    print('##################################')
    
class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone("main.map")

  def initial(self):
    # Creation if initial predator & prey map
    # Prey:
    self.prey = uniform(1)
    # Set probability for prey to 10 % 
    self.prey = self.prey<initialPreyDensity # to be made random .. 
    # Save as map
    self.report(self.prey,"preyMap")
    
    # Predator:
    self.predatormap = uniform(1)
    # Set probability for predator to 5 % 
    self.predator = self.predatormap<0.3
    # Save as map
    self.report(self.predator,"predatorMap")
    
    # infected Predator:
 
    # Set probability for predator to 5 % 
    self.infectedPredator =     pcrand(self.predatormap>initialPredDensity, self.predatormap<initialPredDensity+0.1)
    # Save as map
    self.report(self.infectedPredator,"infectedPredatorMap")

    # Equilibrium state variables:
    self.iter = 0
    self.emptymap = scalar((uniform(1)))*0
    self.preyDensity = [0]*30
    self.predDensity = [0]*30
    self.preyInitDensity = areaaverage(scalar(self.prey),nominal(self.emptymap))
    self.predInitDensity = areaaverage(scalar(self.predator),nominal(self.emptymap))
    self.alreadyPrey = False
    self.alreadyPred = False
    self.preyFinalDensity = 0 
    self.predFinalDensity = 0 
    self.infPredFinalDensity = 0 
    

  def dynamic(self):
    # Determine cells which were occupied by both, prey and predator in last timestep
    both = pcrand(self.prey, self.predator)
        
    # Determine cells which are occupied by prey only in last timestep
    alivePrey = pcrand(self.prey, pcrnot(self.predator))

    # Determine the number of cells with prey only in Von Neumann neighborhood of each cell in last timestep
    preyBirth = window4total(scalar(alivePrey))
    # Determine whether a cell had any cells with prey only in Von Neumann neighborhood or was prey only in the last timestep
    preyBirth = pcror(preyBirth>0.0, alivePrey)
    # All of the cells in preyBirth become prey in the current timestep
    self.prey = preyBirth

    # Determine number of infected predators in neighbourhood
    infectedNeighbour = window4total(scalar(self.infectedPredator))
    infectedNeighbour = infectedNeighbour>0.0
    # Determine cells that do fully reproduce (only have non infected neighbours and can eat)
    doIReproduce = pcrand(both, pcrnot(infectedNeighbour))
    # Determine the number of cells with both, prey and (sound) predator in Von Neumann neighborhood of each cell in last timestep
    birthByNeighbour = window4total(scalar(doIReproduce))
    birthByNeighbour = birthByNeighbour>0.0
    # Determine whether a cell has any cells with both, prey and predator in Von Neumann neighborhood or is itself both prey and predator in last timestep
    predatorBirth = pcror(birthByNeighbour, doIReproduce)
    

    # Determine the cells that become infected  
    infectedCells = pcrand(infectedNeighbour, both)
    self.infectedPredator = infectedCells
    
    # All of the cells in predatorBirth become predator in the current timestep
    self.predator = pcrand(predatorBirth, pcrnot(infectedCells))

    # Save prey and predator distribution of current timestep each as a map
    self.report(self.infectedPredator, 'infPred')
    self.report(self.predator,'predator')
    self.report(self.prey,"prey")
    

    # <--------------------------------- Equilibrium --------------------------------->
    
    # Find equilibrium state. Assumption: equilibrium = mean density of prey/predator for last 10 iterations is in one standard error deviation from the mean density of the previous 30 iterations

    # <------------------------------------ Prey ------------------------------------->

    # Create variables
    preyMean10 = 0
    preyMean30 = 0
    preySeries = 0
    # Save the density of prey for each of the last 30 iterations
    if self.iter<30:
        self.preyDensity[self.iter] = areaaverage(scalar(self.prey),nominal(self.emptymap))
    else:
        for i in range(29):
            self.preyDensity[i]=self.preyDensity[i+1]
        self.preyDensity[29] = areaaverage(scalar(self.prey),nominal(self.emptymap))

    # After the first 30 iterations look for equilibrium
    if self.iter>=29:
        # Calculate the average number of prey in the last 10 iterations 
        for i in range(10):
            preyMean10 = preyMean10 + self.preyDensity[29-i]
        preyMean10 = preyMean10 / 10      

        # Calculate the average number of prey in the last 30 iterations
        for i in range(30):
            preyMean30 = preyMean30 + self.preyDensity[i]
        preyMean30 = preyMean30 / 30

        # Calculate the standard deviation of prey for the last 30 iterations
        for i in range(30):
            preySeries = preySeries + (self.preyDensity[i] - preyMean30)**2
        preySD30 = sqrt(preySeries/30)

        # Calculate the standard deviation of prey for the last 10 iterations
        for i in range(10):
            preySeries = preySeries + (self.preyDensity[29-i] - preyMean10)**2
        preySD10 = sqrt(preySeries/10)

        # Validate equilibrium:
        equlibriumPrey = (preyMean10 > (preyMean30 - preySD30)) & (preyMean10 < (preyMean30 + preySD30))
        
        # Save density of prey in the state of equilibrium
        preyFinalDensity = ifthen(equlibriumPrey, areaaverage(scalar(self.prey),nominal(self.emptymap)))
        self.report(preyFinalDensity,"preyDens")
        if(pcraster.cellvalue(equlibriumPrey,1)[0]):
            # When an equilibrium has already been found in this execution do nothing
            if(not(self.alreadyPrey)):
                    self.preyFinalDensity = pcraster.cellvalue(preyFinalDensity,1)[0]
                    # Set already to true to not calculate any further
                    self.alreadyPrey = True

    # <---------------------------------- Predator ----------------------------------->

    # Create variables
    predMean10 = 0
    predMean30 = 0
    predSeries = 0
    
    # Save the density of pred for each of the last 30 iterations
    if self.iter<30:
        self.predDensity[self.iter] = areaaverage(scalar(self.predator),nominal(self.emptymap))
    else:
        for i in range(29):
            self.predDensity[i]=self.predDensity[i+1]
        self.predDensity[29] = areaaverage(scalar(self.predator),nominal(self.emptymap))

    # After the first 30 iterations look for equilibrium
    if self.iter>=29:
        # Calculate the average number of pred in the last 10 iterations 
        for i in range(10):
            predMean10 = predMean10 + self.predDensity[29-i]
        predMean10 = predMean10 / 10

        # Calculate the average number of pred in the last 30 iterations
        for i in range(30):
            predMean30 = predMean30 + self.predDensity[i]
        predMean30 = predMean30 / 30

        # Calculate the standard deviation of pred for the last 30 iterations
        for i in range(30):
            predSeries = predSeries + (self.predDensity[i] - predMean30)**2
        predSD30 = sqrt(predSeries/30)

        # Calculate the standard deviation of pred for the last 10 iterations
        for i in range(10):
            predSeries = predSeries + (self.predDensity[29-i] - predMean10)**2
        predSD10 = sqrt(predSeries/10)

        # Validate equilibrium:
        equlibriumPred = (predMean10 > (predMean30 - predSD30)) & (predMean10 < (predMean30 + predSD30))
        # Save density of pred in the state of equilibrium
        predFinalDensity = ifthen(equlibriumPred, areaaverage(scalar(self.predator),nominal(self.emptymap)))
        infPredFinalDensity = ifthen(equlibriumPred, areaaverage(scalar(self.infectedPredator),nominal(self.emptymap)))

        predInitDensity = pcraster.cellvalue(self.predInitDensity,1)[0]

        self.report(predFinalDensity,"predDens")
        # Confirm if equilibrium is validated (map has value true)
        if(pcraster.cellvalue(equlibriumPred,1)[0]):
            # When an equilibrium has already been found in this execution do nothing
            if(not(self.alreadyPred)):
                    self.predFinalDensity = pcraster.cellvalue(predFinalDensity,1)[0]
                    self.infPredFinalDensity = pcraster.cellvalue(infPredFinalDensity,1)[0]

                    # Set already to true to not calculate any further
                    self.alreadyPred = True
            
    self.iter=self.iter+1
# Run the model 100 times
nrOfTimeSteps=50
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
i = 0

# preyInitDensity,predInitDensity,self.preyFinalDensityCSV[0],self.predFinalDensityCSV[0]
avg_prey = 0
avg_pred = 0
avg_inf = 0 
while(initialPreyDensity<=1):
    printStatus()
    while(initialPredDensity<=1):
        while(i<int(numberOfTimesExecuted)):
            dynamicModel.run()
            avg_prey += myModel.preyFinalDensity
            avg_pred += myModel.predFinalDensity
            avg_inf += myModel.infPredFinalDensity
            print(i+1)
            print("Averages:")
            print(str(avg_prey/int(i+1))+'/'+str(avg_pred/int(i+1))+"/"+str(avg_inf/int(i+1)))
            i+=1
        avg_prey = avg_prey/int(numberOfTimesExecuted)
        avg_pred = avg_pred/int(numberOfTimesExecuted)
        avg_inf = avg_inf/int(numberOfTimesExecuted)
        print("Average of prey:"+str(avg_prey))
        print("Average of pred:"+str(avg_pred))
        print("Average of infected predator:"+str(avg_inf))
        printToCSV(avg_prey,avg_pred,avg_inf)
        initialPredDensity+=0.1
        print("New distribution values:"+str(initialPreyDensity)+"/"+str(initialPredDensity))
        print()
        avg_prey = 0
        avg_pred = 0
        avg_inf = 0 
        i = 0
    initialPreyDensity+=0.1
    initialPredDensity = 0.0


