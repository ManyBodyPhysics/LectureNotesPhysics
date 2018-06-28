class PairingBasisGen:
    def generateFermi(self):
        for i in range(0,self.nParticles):
            self.below_fermi.append(i)
        for j in range(self.nParticles, self.nSpStates):
            self.above_fermi.append(j)

    def generateStates(self):
        for sp in range(0,self.nSpStates/2):
            self.states.append((sp+1,1))
            self.states.append((sp+1,-1))

    def __init__(self, statesIn, particlesIn):
        self.nSpStates = statesIn
        self.nParticles = particlesIn
        self.below_fermi = []
        self.above_fermi = []
        self.states = []
        self.generateFermi()
        self.generateStates()



