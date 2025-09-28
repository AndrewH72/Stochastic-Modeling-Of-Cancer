generateTumor <- function(sampleSize=100, birthRateS=1, birthRateR=1.2, deathRateS=0.5, deathRateR=0.3, mutationRate1=0.0001, transformationRate=0.0001){
    X = cbind(Time=c(0), S=c(sampleSize), R=c(0))
    updateMatrix = matrix(c(1, 0, -1, 0, 0, 1, -1, 1, 0, 1, 0, -1), nrow=6, ncol=2, byrow=TRUE)
    
    i = 1
    sink(file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/logs/tumorLog.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".transformation_", transformationRate, ".mutation1_", mutationRate1, ".log", sep = ""))
    cat(colnames(X), "\n")
    while(X[i, 3] < sampleSize * 0.1){
        cat(X[i, ], "\n")
        rates = c(birthRateS*X[i, 2], deathRateS*X[i, 2], mutationRate1*X[i, 2], transformationRate*X[i, 2], birthRateR*X[i, 3], deathRateR*X[i, 3])
        r = sum(rates)
        tstep = rexp(1, r)
        
        p = rates/r
        u = runif(1)
        j = which(u <= cumsum(p)) [1]
        X = rbind(X, c(X[i, 1] + tstep, X[i, 2:3] + updateMatrix[j, ]))
        i = i + 1
    }
    cat(tail(X, 1))
    sink()
    
    png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/tumorPlot.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".transformation_", transformationRate, ".mutation1_", mutationRate1, ".png", sep = "")) 
    par(mfrow = c(1, 2))
    plot(X[, 1], X[, 2], ylim=c(0, max(X[, 2])), type="l", lwd=2, main=paste("Susceptible Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab=" Susceptible Population")
    plot(X[, 1], X[, 3], ylim=c(0, max(X[, 3])), type="l", lwd=2, main=paste("Resistant Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population")
    dev.off()
    
    write.csv(X, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/tumorData.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".transformation_", transformationRate, ".mutation1_", mutationRate1, ".csv", sep = ""), row.names=FALSE)
    return(X)
}

cyclicTreatment <- function(sampleSize=100, birthRateS=1, birthRateR=0.8, deathRateS=0.5, deathRateR=0.3, mutationRate1=0.0001, transformationRate=0.0001, drugInducedDeathRate=0.6, cycleDuration=2.5, maxDuration=5){
    X = tail(generateTumor(sampleSize, birthRateS, birthRateR, deathRateS, birthRateR, mutationRate1, transformationRate), 1)
    colnames(X) = c("Time", "S", "R")
    updateMatrix = matrix(c(1, 0, -1, 0, 0, 1, -1, 1, 0, 1, 0, -1), nrow=6, ncol=2, byrow=TRUE)
    i = 1
    
    sink(file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/logs/treatmentLog.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".drugDeath_", drugInducedDeathRate,".transformation_", transformationRate, ".mutation1_", mutationRate1, ".log", sep=""))
    cat(colnames(X), "\n")
    while(X[i, 1] < X[1, 1] + maxDuration & (X[i, 2] + X[i, 3]) > 0){
        cat(X[i, ], "\n")
        rates = c(birthRateS*X[i, 2], (deathRateS + drugInducedDeathRate)*X[i, 2], mutationRate1*X[i, 2], transformationRate*X[i, 2], birthRateR*X[i, 3], (deathRateR + drugInducedDeathRate)*X[i, 3])
        r = sum(rates)
        tstep = rexp(1, r)
        
        p = rates/r
        u = runif(1)
        j = which(u <= cumsum(p)) [1]
        X = rbind(X, c(X[i, 1] + tstep, X[i, 2:3] + updateMatrix[j, ]))
        i = i + 1
    }
    cat(tail(X, 1))
    sink()
    
    png(filename=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/treatmentPlot.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".drugDeath_", drugInducedDeathRate,".transformation_", transformationRate, ".mutation1_", mutationRate1, ".png", sep="")) 
    par(mfrow = c(1, 2))
    plot(X[, 1], X[, 2], ylim=c(0, max(X[, 2])), type="l", lwd=2, main=paste("Susceptible Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Susceptible Population")
    plot(X[, 1], X[, 3], ylim=c(0, max(X[, 3])), type="l", lwd=2, main=paste("Resistant Population vs Time for Starting Size", sampleSize), cex.main=0.70, xlab="Time", ylab="Resistant Population")
    dev.off()
    
    write.csv(X, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/treatmentData.sampleSize_", sampleSize, ".birthS_", birthRateS, ".birthR_", birthRateR, ".deathS_", deathRateS, ".deathR_", deathRateR, ".drugDeath_", drugInducedDeathRate,".transformation_", transformationRate, ".mutation1_", mutationRate1, ".csv", sep=""), row.names=FALSE)
    return(X)
}

cyclicTreatment()
