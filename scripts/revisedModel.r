################
# Komarova Model
################
komarovaModel2Drugs = function(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, maxTime){
    ##############################################################
    # Data frame that contains the different cell population sizes.
    ##############################################################
    tumorData = data.frame(
        Time = 0,
        S = sampleSize, # Susceptible 
        R1 = 0, # Resistant to Drug 1
        R2 = 0, # Resistant to Drug 2
        RB = 0 # Resistant to Drug 1 and Drug 2
    )
    
    #######################################################################
    ## Different possible events and the ways cell populations get updated. 
    #######################################################################
    updateMatrix = matrix(c(
        # Faithful Reproduction
        1, 0, 0, 0,
        
        # Mutation
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
        
        # Death
        -1, 0, 0, 0,
        0, -1, 0, 0,
        0, 0, -1, 0,
        0, 0, 0, -1,
        
        # Transformation S -> R
        -1, 1, 0, 0,
        -1, 0, 1, 0,
        
        # Transformation R -> R
        0, -1, 0, 1,
        0, 0, -1, 1
        
    ), ncol=4, byrow=TRUE)
    colnames(updateMatrix) = c("S", "R1", "R2", "RB")
    rownames(updateMatrix) = c("Faithful Division", "Mutation R1", "Mutation R2", "Mutation RB", "Death S", "Death R1", "Death R2", "Death RB", "Transformation S R1", "Transformation S R2", "Transformation R1 RB", "Transformation R2 RB")
    
    #####################
    # Important Variables
    #####################
    iterations = 0
    time = tumorData$Time
    currentState = c(S=tumorData$S, R1=tumorData$R1, R2=tumorData$R2, RB=tumorData$RB)
    sumOfMutations = mutationRate1 + mutationRate2 + mutationRateB
    sumOfDrugInducedDeathRates = drugInducedDeathRate1 + drugInducedDeathRate2
    
    
    ###########
    # Main Loop
    ###########
    sink(file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/logs/RevisedModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".log", sep = ""))
    cat(colnames(tumorData), "\n")
    while(time < maxTime & currentState["RB"] < 1){
        cat(as.character(tumorData[iterations + 1, ]), "\n")
        currentState[currentState < 0] = 0
        
        rates = c(
        # Faithful Reproduction
        birthRate * (1 - sumOfMutations) * currentState["S"],
    
        # Mutation
        birthRate * (1 - mutationRateB) * currentState["R1"],
        birthRate * (1 - mutationRateB) * currentState["R2"],
        birthRate * currentState["RB"],
    
        # Death
        (deathRate + sumOfDrugInducedDeathRates) * currentState["S"],
        (deathRate + drugInducedDeathRate2) * currentState["R1"],
        (deathRate + drugInducedDeathRate1) * currentState["R2"],
        deathRate * currentState["RB"],
    
        # Transformation S -> R
        birthRate * mutationRate1 * currentState["S"],
        birthRate * mutationRate2 * currentState["S"],
    
        # Transformation R -> R
        birthRate * mutationRateB * currentState["R1"],
        birthRate * mutationRateB * currentState["R2"]
        )
        
        totalRates = sum(rates)
        if(totalRates == 0){
            break
        }
        
        # Identify the event that will take place at the time step. 
        tstep = rexp(1, totalRates)
        event = runif(1, 0, totalRates)
        eventIdx = which.max(event <= cumsum(rates))
        
        # Update the state.
        currentState = currentState + updateMatrix[eventIdx, ]
        newRow = data.frame(
            Time = time + tstep,
            S = currentState["S"],
            R1 = currentState["R1"],
            R2 = currentState["R2"],
            RB = currentState["RB"]
        )
        tumorData = rbind(tumorData, newRow)
        rownames(tumorData) = NULL
        time = time + tstep
        iterations = iterations + 1
    }
    cat(as.character(tail(tumorData, 1)))
    sink()
    
    png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/RevisedModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".png", sep = "")) 
    par(mfrow = c(2, 2))
    plot(tumorData[, 1], tumorData[, 2], ylim=c(0, max(tumorData[, 2])), type="l", lwd=2, main=paste("Susceptible Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab=" Susceptible Population")
    plot(tumorData[, 1], tumorData[, 3], ylim=c(0, max(tumorData[, 3])), type="l", lwd=2, main=paste("Resistant Population 1 vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population 1")
    plot(tumorData[, 1], tumorData[, 4], ylim=c(0, max(tumorData[, 4])), type="l", lwd=2, main=paste("Resistant Population 2 vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population 2")
    plot(tumorData[, 1], tumorData[, 5], ylim=c(0, max(tumorData[, 5])), type="l", lwd=2, main=paste("Resistant Population B vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population B")
    dev.off()
    
    write.csv(tumorData, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/RevisedModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".csv", sep = ""), row.names=FALSE)
    cat("Max Susceptible Population:", max(tumorData$S), "\n", "Max Resistant1 Population:", max(tumorData$R1), "\n", "Max Resistant2 Population:", max(tumorData$R2), "\n", "Max ResistantB Population:", max(tumorData$RB), "\n")
    return(tail(tumorData$Time, n=1))
}


##############
# Cyclic Model
##############
cyclicModel = function(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, cycleTime1, cycleTime2, totalCycles){
    tumorData = data.frame(
        Time = 0,
        S = sampleSize, # Susceptible 
        R1 = 0, # Resistant to Drug 1
        R2 = 0, # Resistant to Drug 2
        RB = 0 # Resistant to Drug 1 and Drug 2
    )
    
    updateMatrix = matrix(c(
        # Faithful Reproduction
        1, 0, 0, 0,
        
        # Mutation
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
        
        # Death
        -1, 0, 0, 0,
        0, -1, 0, 0,
        0, 0, -1, 0,
        0, 0, 0, -1,
        
        # Transformation S -> R
        -1, 1, 0, 0,
        -1, 0, 1, 0,
        
        # Transformation R -> R
        0, -1, 0, 1,
        0, 0, -1, 1
        
    ), ncol=4, byrow=TRUE)
    colnames(updateMatrix) = c("S", "R1", "R2", "RB")
    rownames(updateMatrix) = c("Faithful Division", "Mutation R1", "Mutation R2", "Mutation RB", "Death S", "Death R1", "Death R2", "Death RB", "Transformation S R1", "Transformation S R2", "Transformation R1 RB", "Transformation R2 RB")
    
    iterations = 0
    time = tumorData$Time
    currentState = c(S=tumorData$S, R1=tumorData$R1, R2=tumorData$R2, RB=tumorData$RB)
    sumOfMutations = mutationRate1 + mutationRate2 + mutationRateB
    sumOfDrugInducedDeathRates = drugInducedDeathRate1 + drugInducedDeathRate2
    
    
    sink(file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/logs/CyclicModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".log", sep = ""))
    cat(colnames(tumorData), "\n")
    for(cycle in 1:totalCycles){
        cycleEnd = time + cycleTime1
        while(time < cycleEnd){
            cat(as.character(tumorData[iterations + 1, ]), "\n")
            currentState[currentState < 0] = 0
            
            rates = c(
                # Faithful Reproduction
                birthRate * (1 - sumOfMutations) * currentState["S"],
                
                # Mutation
                birthRate * (1 - mutationRateB) * currentState["R1"],
                birthRate * (1 - mutationRateB) * currentState["R2"],
                birthRate * currentState["RB"],
                
                # Death
                (deathRate + sumOfDrugInducedDeathRates) * currentState["S"],
                deathRate * currentState["R1"],
                (deathRate + drugInducedDeathRate1) * currentState["R2"],
                deathRate * currentState["RB"],
                
                # Transformation S -> R
                birthRate * mutationRate1 * currentState["S"],
                birthRate * mutationRate2 * currentState["S"],
                
                # Transformation R -> R
                birthRate * mutationRateB * currentState["R1"],
                birthRate * mutationRateB * currentState["R2"]
            )
            
            totalRates = sum(rates)
            if(totalRates == 0){
                break
            }
            
            # Identify the event that will take place at the time step. 
            tstep = rexp(1, totalRates)
            event = runif(1, 0, totalRates)
            eventIdx = which.max(event <= cumsum(rates))
            
            # Update the state.
            currentState = currentState + updateMatrix[eventIdx, ]
            newRow = data.frame(
                Time = time + tstep,
                S = currentState["S"],
                R1 = currentState["R1"],
                R2 = currentState["R2"],
                RB = currentState["RB"]
            )
            tumorData = rbind(tumorData, newRow)
            rownames(tumorData) = NULL
            time = time + tstep
            iterations = iterations + 1
        }
        
        cycleEnd = time + cycleTime2
        while(time < cycleEnd){
            cat(as.character(tumorData[iterations + 1, ]), "\n")
            currentState[currentState < 0] = 0
            
            rates = c(
                # Faithful Reproduction
                birthRate * (1 - sumOfMutations) * currentState["S"],
                
                # Mutation
                birthRate * (1 - mutationRateB) * currentState["R1"],
                birthRate * (1 - mutationRateB) * currentState["R2"],
                birthRate * currentState["RB"],
                
                # Death
                (deathRate + sumOfDrugInducedDeathRates) * currentState["S"],
                (deathRate + drugInducedDeathRate2) * currentState["R1"],
                deathRate * currentState["R2"],
                deathRate * currentState["RB"],
                
                # Transformation S -> R
                birthRate * mutationRate1 * currentState["S"],
                birthRate * mutationRate2 * currentState["S"],
                
                # Transformation R -> R
                birthRate * mutationRateB * currentState["R1"],
                birthRate * mutationRateB * currentState["R2"]
            )
            
            totalRates = sum(rates)
            if(totalRates == 0){
                break
            }
            
            # Identify the event that will take place at the time step. 
            tstep = rexp(1, totalRates)
            event = runif(1, 0, totalRates)
            eventIdx = which.max(event <= cumsum(rates))
            
            # Update the state.
            currentState = currentState + updateMatrix[eventIdx, ]
            newRow = data.frame(
                Time = time + tstep,
                S = currentState["S"],
                R1 = currentState["R1"],
                R2 = currentState["R2"],
                RB = currentState["RB"]
            )
            tumorData = rbind(tumorData, newRow)
            rownames(tumorData) = NULL
            time = time + tstep
            iterations = iterations + 1
        }
    }
    cat(as.character(tail(tumorData, 1)))
    sink()
    
    png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/CyclicModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".png", sep = "")) 
    par(mfrow = c(2, 2))
    plot(tumorData[, 1], tumorData[, 2], ylim=c(0, max(tumorData[, 2])), type="l", lwd=2, main=paste("Susceptible Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab=" Susceptible Population")
    plot(tumorData[, 1], tumorData[, 3], ylim=c(0, max(tumorData[, 3])), type="l", lwd=2, main=paste("Resistant Population 1 vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population 1")
    plot(tumorData[, 1], tumorData[, 4], ylim=c(0, max(tumorData[, 4])), type="l", lwd=2, main=paste("Resistant Population 2 vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population 2")
    plot(tumorData[, 1], tumorData[, 5], ylim=c(0, max(tumorData[, 5])), type="l", lwd=2, main=paste("Resistant Population B vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population B")
    dev.off()
    
    write.csv(tumorData, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/CyclicModel/tumorLog.sampleSize_", sampleSize, ".birthRate_", birthRate, ".deathRate_", deathRate, ".mutationRate1_", mutationRate1, ".mutationRate2_", mutationRate2, ".mutationRateRB_", mutationRateB, ".drugInducedDeathRate1_", drugInducedDeathRate1, ".drugInducedDeathRate2_", drugInducedDeathRate2, ".dateTime_", Sys.time(), ".csv", sep = ""), row.names=FALSE)
    
}


#############
# Main Script
#############
mainScript = function(){
    maxTime = 100
    maxSimulations = 10
    # Birth and Death are arbitrarily chosen, all that matters is the turnover rate (Death/Birth) is between 0.1-0.5
    birthRate = 0.99
    deathRate = 0.5
    # birthRate = c(1.5, 1.1, 1, 0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
    # deathRate = c(0.75, 0.55, 0.5, 0.5, 0.475, 0.45, 0.42, 0.4, 0.375, 0.35)
    
    # Mutation rates were arbitrarily chosen. Komarova provided a range of mutation rates to try out.
    mutationRate1 = 10e-4 
    mutationRate2 = 10e-4
    mutationRateB = 10e-4
    
    # Komarova states there are many ways to simulate the drug induced death rate, in her model, she simply sets each drug's induced death rate to the treatment intensity if cells are susceptible to at least one drug and 0 if they are susceptible to 0.
    # Formula for treatment intensity was taken from Section 4.1 The case of one drug.
    treatmentIntensity = 2 * (birthRate - deathRate)
    drugInducedDeathRate1 = 0.01
    drugInducedDeathRate2 = 0.01
    
    # Initial population was arbitrarily chosen. It appears Komarova used 100 for both of her numerical simulations, however, these were simply used to test the dependence between variables, so it may not be accurate.
    sampleSize = 1000
    
    avgTimeVec = c()
    for(i in 1:maxSimulations){
        avgTime = komarovaModel2Drugs(birthRate, deathRate, mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, maxTime)
        # avgTime = komarovaModel2Drugs(birthRate[i], deathRate[i], mutationRate1, mutationRate2, mutationRateB, drugInducedDeathRate1, drugInducedDeathRate2, sampleSize, maxTime)
        avgTimeVec = c(avgTimeVec, avgTime)
    }
    print(mean(avgTimeVec))
}

#################
# Run Main Script
################
mainScript()
