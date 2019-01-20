from pcraster import *
from pcraster.framework import *

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone("main.map")

  def initial(self):
      # Creation if initial predator & prey map 
    self.prey = uniform(1)
    self.prey = self.prey<0.11
    self.report(self.prey,"preyMap")
    self.predator = uniform(1)
    self.predator = self.predator<0.11
    self.report(self.predator,"predatorMap")

    print("Finished building maps for prey&predator distribution")
  def dynamic(self):
    # Determine cells which are occupied by prey and predator
    preyFeeds = pcrand(self.prey,self.predator)
    self.report(preyFeeds,"fed")

    # Determine if a mother for offspring is available in next cell
    preyHasMom = window4total(scalar(self.prey))
    preyHasMom = preyHasMom>0.0
    self.report(preyHasMom,"mom")
    predatorHasMom = window4total(scalar(self.predator))
    predatorHasMom = predatorHasMom>0.0
    self.report(predatorHasMom,'mompred')
    # Combine both results to determine if a cell is occupied by prey 
    prey = pcrand(preyHasMom,pcrnot(preyFeeds))

    predator = pcror(pcrand(self.prey,self.predator),predatorHasMom)
    self.predator = predator
    self.report(predator,'predator')
    self.prey = prey
    self.report(prey,"prey")

nrOfTimeSteps=10
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()

  




