generateTumor <- function(sampleSize=100, birthRate=1, deathRate=0.5, mutationRate1=0.0001, transformationRate=0.0001){
    X = cbind(Time=c(0), S=c(sampleSize), R=c(0))
    updateMatrix = matrix(c(1, 0, -1, 0, 0, 1, -1, 1), nrow=4, ncol=2, byrow=TRUE)
    
    i = 1
    while(X[i, 3] < sampleSize * 0.1){
        rates = c(birthRate*X[i, 2], deathRate*X[i, 2], mutationRate1*X[i, 2], transformationRate*X[i, 2])
        r = sum(rates)
        tstep = rexp(1, r)
        
        p = rates/r
        u = runif(1)
        j = which(u <= cumsum(p)) [1]
        X = rbind(X, c(X[i, 1] + tstep, X[i, 2:3] + updateMatrix[j, ]))
        i = i + 1
    }
    write.csv(X, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/cellData", sampleSize, ".birth_", birthRate, ".death_", deathRate, ".transformation_", transformationRate, ".mutation1_", mutationRate1, ".csv", sep=""), row.names=FALSE)
    return(X)
}

cyclicTreatment <- function(preTreatmentTumor, birthRate=1, deathRate=0.5, mutationRate1=0.0001, transformationRate=0.0001, drugInducedDeathRate=0.25, cycleDuration=12, maxDuration=60){
    X = matrix(preTreatmentTumor, nrow=1)
    colnames(X) = c("Time", "S", "R")
    i = 1
    
    while(X[i, 1] < X[1, 1] + maxDuration & (X[i, 2] + X[i, 3]) > 0){
        rates = c(birthRate*X[i, 2], (deathRate + drugInducedDeathRate)*X[i, 2], mutationRate1*X[i, 2], transformationRate*X[i, 2])
        r = sum(rates)
        tstep = rexp(1, r)
        
        p = rates/r
        u = runif(1)
        j = which(u <= cumsum(p)) [1]
        X = rbind(X, c(X[i, 1] + tstep, X[i, 2:3] + updateMatrix[j, ]))
        i = i + 1
    }
    write.csv(X, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/cellData", sampleSize, ".birth_", birthRate, ".death_", deathRate, ".drugDeath_", drugInducedDeathRate,".transformation_", transformationRate, ".mutation1_", mutationRate1, ".csv", sep=""), row.names=FALSE)
    return(X)
}

preTreatmentTumor = tail(generateTumor(sampleSize=100), 1)
