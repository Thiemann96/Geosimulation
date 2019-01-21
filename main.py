from pcraster import *
from pcraster.framework import *

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone("main.map")

  def initial(self):
    # Creation if initial predator & prey map 
    self.prey = uniform(1)
    self.prey = self.prey<0.10
    self.report(self.prey,"preyMap")
    self.predator = uniform(1)
    self.predator = self.predator<0.10
    self.report(self.predator,"predatorMap")

    print("Finished building maps for prey&predator distribution")
  def dynamic(self):
    # Determine cells which are occupied by prey and predator
    both = pcrand(self.prey,self.predator)
    self.report(both,"pboth")
    alivePrey = pcrand(self.prey, pcrnot(self.predator))
    self.report(alivePrey,"aPrey")

    # Determine if a mother for offspring is available in next cell
    preyBirthRate = window4total(scalar(alivePrey))
    preyBirthRate = preyBirthRate>0.0
    preyBirth = pcror(preyBirthRate, alivePrey)
    self.report(preyBirth,"pyBirth")
    predatorBirth = window4total(scalar(both))
    predatorBirth = pcror(predatorBirth>0.0, both)
    self.report(predatorBirth,'prBirth')
    # Combine both results to determine if a cell is occupied by prey 
    prey = preyBirth

    predator = predatorBirth
    self.predator = predator
    self.report(predator,'predator')
    self.prey = prey
    self.report(prey,"prey")

nrOfTimeSteps=100
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()


