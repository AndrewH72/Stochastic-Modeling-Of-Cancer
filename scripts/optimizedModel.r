optimizedKomarova2D = function(birthRate=1.0, deathRate=0.5, mutationRate1=10e-4, mutationRate2=10e-4, mutationRateB=10e-4, drugInducedDeathRate1=0.1, drugInducedDeathRate2=0.1, sampleSize=1000, maxTime=100, tau=0.1){
   estimatedSteps = (maxTime / tau) * 2
   tumorHistory = matrix(0, nrow=estimatedSteps, ncol=5)
   colnames(tumorHistory) = c("Time", "S", "R1", "R2", "RB")
   
   time = 0   
   iter = 1
   currentState = c(S=sampleSize, R1=0, R2=0, RB=0)
   tumorHistory[iter, ] = c(time, currentState)
   
   sumOfMutations = mutationRate1 + mutationRate2 + mutationRateB
   sumOfDrugInducedDeathRates = drugInducedDeathRate1 + drugInducedDeathRate2
   
   updateMatrix = matrix(c(
       1, 0, 0, 0, # Faithful S
       0, 1, 0, 0, # Mutation R1
       0, 0, 1, 0, # Mutation R2
       0, 0, 0, 1, # Mutation RB
       -1, 0, 0, 0, # Death S
       0, -1, 0, 0, # Death R1
       0, 0, -1, 0, # Death R2
       0, 0, 0, -1, # Death RB
       -1, 1, 0, 0, # Trans S -> R1
       -1, 0, 1, 0, # Trans S -> R2
       0, -1, 0, 1, # Trans R1 -> RB
       0, 0, -1, 1 # Trans R2 -> RB
   ), ncol=4, byrow=TRUE)
   
   while(time < maxTime & currentState["RB"] < 1){
       
   }
   
}