from pcraster import *
from pcraster.framework import *

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone("main.map")

  def initial(self):
    # Creation if initial predator & prey map
    # Prey:
    self.prey = uniform(1)
    # Set probability for prey to 10 % 
    self.prey = self.prey<0.10
    # Save as map
    self.report(self.prey,"preyMap")
    
    # Predator:
    self.predator = uniform(1)
    # Set probability for predator to 10 % 
    self.predator = self.predator<0.10
    # Save as map
    self.report(self.predator,"predatorMap")

    print("Finished building maps for prey&predator distribution")
  def dynamic(self):
    # Determine cells which were occupied by both, prey and predator in last timestep
    both = pcrand(self.prey,self.predator)
        
    # Determine cells which are occupied by prey only in last timestep
    alivePrey = pcrand(self.prey, pcrnot(self.predator))

    # Determine the number of cells with prey only in Von Neumann neighborhood of each cell in last timestep
    preyBirth = window4total(scalar(alivePrey))
    # Determine whether a cell had any cells with prey only in Von Neumann neighborhood or was prey only in the last timestep
    preyBirth = pcror(preyBirth>0.0, alivePrey)
    # All of the cells in preyBirth become prey in the current timestep
    self.prey = preyBirth
   
    # Determine the number of cells with both, prey and predator in Von Neumann neighborhood of each cell in last timestep
    predatorBirth = window4total(scalar(both))
    # Determine whether a cell has any cells with both, prey and predator in Von Neumann neighborhood or is itself both prey and predator in last timestep
    predatorBirth = pcror(predatorBirth>0.0, both)
    # All of the cells in predatorBirth become predator in the current timestep
    self.predator = predatorBirth
    
    # Save prey and predator distribution of current timestep each as a map
    self.report(predator,'predator')
    self.report(prey,"prey")
    
# Run the model 100 times
nrOfTimeSteps=100
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()


