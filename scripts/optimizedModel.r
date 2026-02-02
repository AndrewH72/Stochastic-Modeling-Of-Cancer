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
           birthRate * mutationRateB * currentState["R2"]
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

optimizedCyclicalTreatment = function(birthRate=1.0, deathRate=0.5, mutationRate1=10e-4, mutationRate2=10e-4, mutationRateB=10e-4, drugInducedDeathRate1=0.1, drugInducedDeathRate2=0.1, sampleSize=1000, tau=0.1, cycleTime1=5, cycleTime2=10, totalCycles=100){
    estimatedSteps = (totalCycles / tau) * 2
    tumorHistory = matrix(0, nrow=estimatedSteps, ncol=5)
    colnames(tumorHistory) = c("Time", "S", "R1", "R2", "RB")
    
    time = 0   
    iter = 1
    currentState = c(S=sampleSize, R1=0, R2=0, RB=0)
    tumorHistory[iter, ] = c(time, currentState)
    
    sumOfMutations = mutationRate1 + mutationRate2 + mutationRateB
    
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
    
    for(cycle in 1:totalCycles){
        while(time < cycleTime1){
            rates = c(
                birthRate * (1 - sumOfMutations) * currentState["S"],
                birthRate * (1 - mutationRateB) * currentState["R1"],
                birthRate * (1 - mutationRateB) * currentState["R2"],
                birthRate * currentState["RB"],
                (deathRate + drugInducedDeathRate1) * currentState["S"],
                deathRate * currentState["R1"],
                (deathRate + drugInducedDeathRate1) * currentState["R2"],
                deathRate * currentState["RB"],
                birthRate * mutationRate1 * currentState["S"],
                birthRate * mutationRate2 * currentState["S"],
                birthRate * mutationRateB * currentState["R1"],
                birthRate * mutationRateB * currentState["R2"]
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
        while(time < cycleTime2){
            rates = c(
                birthRate * (1 - sumOfMutations) * currentState["S"],
                birthRate * (1 - mutationRateB) * currentState["R1"],
                birthRate * (1 - mutationRateB) * currentState["R2"],
                birthRate * currentState["RB"],
                (deathRate + drugInducedDeathRate2) * currentState["S"],
                (deathRate + drugInducedDeathRate2) * currentState["R1"],
                deathRate * currentState["R2"],
                deathRate * currentState["RB"],
                birthRate * mutationRate1 * currentState["S"],
                birthRate * mutationRate2 * currentState["S"],
                birthRate * mutationRateB * currentState["R1"],
                birthRate * mutationRateB * currentState["R2"]
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
    }
    tumorHistory = tumorHistory[1:iter, ]
    tumorData = as.data.frame(tumorHistory)
    return(tumorData)
}

runSimulationKomarova = function(){
    set.seed(42)
    maxTime = 10000
    deathRates = 0.99
    # deathRates = seq(0.5, 0.9, 0.05)
    turnOverRatios = c(0.99)
    # turnOverRatios = seq(0.5, 0.9, 0.05)
    mutationRates1 = c(10e-6)
    mutationRates2 = c(10e-6)
    mutationRatesB = c(10e-6)
    drugInducedDeathRates1 = c(0.05)
    drugInducedDeathRates2 = c(0.05)
    # mutationRates1 = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    # mutationRates2 = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    # mutationRatesB = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    sampleSize = 1000
    # drugInducedDeathRates1 = seq(0.01, 0.1, 0.01)
    # drugInducedDeathRates2 = seq(0.01, 0.1, 0.01)
    # tau = 0.1
    taus = 10^(-4:-1)
    
    for(turnOverRatio in turnOverRatios){
        for(deathRate in deathRates){
            birthRate = deathRate / turnOverRatio
            for(mutationRate1 in mutationRates1){
                for(mutationRate2 in mutationRates2){
                    for(mutationRateB in mutationRatesB){
                        for(drugInducedDeathRate1 in drugInducedDeathRates1){
                            for(drugInducedDeathRate2 in drugInducedDeathRates2){
                                for(tau in taus){
                                    fileName = paste("birthRate.", birthRate, "_deathRate.", deathRate, "_mutationRate1.", mutationRate1, "_mutationRate2.", mutationRate2, "_mutationRateB.", mutationRateB, "_drugInducedDeathRate1.", drugInducedDeathRate1, "_drugInducedDeathRate2.", drugInducedDeathRate2, "_sampleSize.", sampleSize, "_tau.", tau, sep="")
                                    startTime = Sys.time()
                                    result = optimizedKomarova2D(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, maxTime, tau)
                                    endTime = Sys.time()
                                    print(paste("Run Time:", endTime - startTime))
                                    
                                    png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/OptimizedModel/tau/", fileName, ".png", sep=""))
                                    par(mfrow=c(2,2))
                                    plot(result$Time, result$S, type="l", main="S", xlab="Time", ylab="Population Size")
                                    plot(result$Time, result$R1, type="l", main="R1", xlab="Time", ylab="Population Size")
                                    plot(result$Time, result$R2, type="l", main="R2", xlab="Time", ylab="Population Size")
                                    plot(result$Time, result$RB, type="l", main="RB", xlab="Time", ylab="Population Size")
                                    dev.off() 
                                    
                                    write.csv(result, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/OptimizedModel/tau/", fileName, ".csv", sep=""), row.names=FALSE)
                                    
                                    sink(file="/Users/andrewhsu/Projects/PREP-NURA/NURA/data/OptimizedModel/statistics.txt", append=TRUE)
                                    print(paste("tau/", fileName, sep=""))
                                    print("Average of Each Column:") 
                                    print(colMeans(result))
                                    
                                    print("Max of Each Column:") 
                                    print(apply(result, 2, max))
                                    cat("\n")
                                    sink()
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


runSimulationKomarova()


runSimulationCyclic = function(){
    totalCycles = 100
    cycleTime1 = 5
    cycleTime2 = 10
    deathRates = c(0.99)
    # deathRates = seq(0.5, 0.9, 0.05)
    turnOverRatios = c(0.99)
    # turnOverRatios = seq(0.5, 0.9, 0.05)
    mutationRates1 = c(10e-6)
    mutationRates2 = c(10e-6)
    mutationRatesB = c(10e-6)
    drugInducedDeathRates1 = c(0.05)
    # drugInducedDeathRates2 = c(0.05)
    # mutationRates1 = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    # mutationRates2 = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    # mutationRatesB = c(10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4)
    sampleSize = 1000
    # drugInducedDeathRates1 = seq(0.01, 0.1, 0.01)
    # drugInducedDeathRates2 = seq(0.01, 0.1, 0.01)
    # taus = c(0.1)
    taus = 10^(-4:-1)
    
    for(turnOverRatio in turnOverRatios){
        for(deathRate in deathRates){
            birthRate = deathRate / turnOverRatio
            for(mutationRate1 in mutationRates1){
                for(mutationRate2 in mutationRates2){
                    for(mutationRateB in mutationRatesB){
                        for(drugInducedDeathRate1 in drugInducedDeathRates1){
                            drugInducedDeathRate2 = drugInducedDeathRate1 / 2
                            for(tau in taus){
                                fileName = paste("birthRate.", birthRate, "_deathRate.", deathRate, "_mutationRate1.", mutationRate1, "_mutationRate2.", mutationRate2, "_mutationRateB.", mutationRateB, "_drugInducedDeathRate1.", drugInducedDeathRate1, "_drugInducedDeathRate2.", drugInducedDeathRate2, "_sampleSize.", sampleSize, "_tau.", tau, "_totalCycles.", totalCycles, "_cycleTime1.", cycleTime1, "_cycleTime2.", cycleTime2, sep="")
                                startTime = Sys.time()
                                result = optimizedCyclicalTreatment(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, tau, cycleTime1, cycleTime2, totalCycles)
                                endTime = Sys.time()
                                print(paste("Run Time:", endTime - startTime))
                                
                                png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/OptimizedCyclic-SSWL/tau/", fileName, ".png", sep=""))
                                par(mfrow=c(2,2))
                                plot(result$Time, result$S, type="l", main="S", xlab="Time", ylab="Population Size")
                                plot(result$Time, result$R1, type="l", main="R1", xlab="Time", ylab="Population Size")
                                plot(result$Time, result$R2, type="l", main="R2", xlab="Time", ylab="Population Size")
                                plot(result$Time, result$RB, type="l", main="RB", xlab="Time", ylab="Population Size")
                                dev.off() 
                                
                                write.csv(result, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/OptimizedCyclic-SSWL/tau/", fileName, ".csv", sep=""), row.names=FALSE)
                                
                                sink(file="/Users/andrewhsu/Projects/PREP-NURA/NURA/data/OptimizedCyclic-SSWL/statistics.txt", append=TRUE)
                                print(paste("SSWL-tau/", fileName, sep=""))
                                print("Average of Each Column:") 
                                print(colMeans(result))
                                
                                print("Max of Each Column:") 
                                print(apply(result, 2, max))
                                cat("\n")
                                sink()
                            }
                        }
                    }
                }
            }
        }
    }
}

runSimulationCyclic()
