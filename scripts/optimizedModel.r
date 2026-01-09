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
       rates = c(
           birthRate * (1 - sumOfMutations) * currentState["S"],
           birthRate * (1 - mutationRateB) * currentState["R1"],
           birthRate * (1 - mutationRateB) * currentState["R2"],
           birthRate * currentState["RB"],
           (deathRate + sumOfDrugInducedDeathRates) * currentState["S"],
           (deathRate + drugInducedDeathRate2) * currentState["R1"],
           (deathRate + drugInducedDeathRate1) * currentState["R2"],
           deathRate * currentState["RB"],
           birthRate * mutationRate1 * currentState["S"],
           birthRate * mutationRate2 * currentState["S"],
           birthRate * mutationRateB * currentState["R1"],
           birthRate * mutationRateB * currentState["R2"],
       )
       
       totalRate = sum(rates)
       if(totalRate == 0) break
       
       currentTau = tau
       numEvents = rpois(length(rates), lambda=rates*currentTau)
       netChange = colSums(updateMatrix * numEvents)
       candidateState = currentState + netChange
       
       if(any(candidateState < 0)){
           candidateState[candidateState < 0] = 0
       }
       
       currentState = candidateState
       time = time + currentTau
       iter = iter + 1
       
       if(iter > nrow(tumorHistory)){
           tumorHistory = rbind(tumorHistory, matrix(0, nrow=estimatedSteps, ncol=5))
       }
       tumorHistory[iter, ] = c(time, currentState)
   }
   tumorHistory = tumorHistory[1:iter, ]
   tumorData = as.data.frame(tumorHistory)
   return(tumorData)
}

runSimulation = function(){
    maxTime = 100
    birthRate = 0.99
    deathRate = 0.5
    mutationRate1 = 10e-4
    mutationRate2 = 10e-4
    mutationRateB = 10e-4
    sampleSize = 1000
    drugInducedDeathRate1 = 0.01
    drugInducedDeathRate2 = 0.01
    tau = 0.1
    
    startTime = Sys.time()
    result = optimizedKomarova2D(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, maxTime, tau)
    endTime = Sys.time()
    
    print(paste("Run Time:", endTime - startTime))
    par(mfrow=c(2,2))
    plot(result$Time, result$S, type="l", main="S", xlab="Time")
    plot(result$Time, result$R1, type="l", main="R1", xlab="Time")
    plot(result$Time, result$R2, type="l", main="R2", xlab="Time")
    plot(result$Time, result$RB, type="l", main="RB", xlab="Time")
}