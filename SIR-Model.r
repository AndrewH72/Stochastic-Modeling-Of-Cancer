####################

## TO DO:
## Turnover rate is currently 0.5 (l=1, d=0.5, h=0.01, u1=0.01)
## Test different l and d rates, closer to 0.1, but keep the turnover rate around 0.1-0.5
## Add in a drug-induced death rate H_s (cyclical treatment)
## One drug for a long time vs Another drug for a long time vs Both drugs alternating

####################

set.seed(123)
n = c(10, 100, 1000)

for(sampleSize in n){
    X = cbind(Time=c(0), S=c(sampleSize), R=c(0)) # store output
    l = 1 # birth rate
    d = 0.5 # death rate
    h = 0.01 # transformation rate
    u1 = 0.01 # mutation rate of drug 1
    
    # Update matrix for each event: faithful division, death, mutation, transformation.  
    updateMatrix = matrix(c(1, 0, -1, 0, 0, 1, -1, 1), nrow=4, ncol=2, byrow=TRUE)
    
    i = 1 # row counter for matrix
    # stochastic simulation algorithm
    while(X[i, 3] < sampleSize * 0.1){
        rates = c(l*X[i, 2], d*X[i, 2], u1*X[i, 2], h*X[i, 2]) # indiv. event rates
        r = sum(rates) # total rate of some event happening
        tstep = rexp(1, r) # time to the next event
        
        p = rates/r # probability vector for each event type
        u = runif(1) # sample from standard uniform
        j = which(u <= cumsum(p)) [1] # sample which event occurred
        
        X = rbind(X, c(X[i, 1] + tstep, X[i, 2:3] + updateMatrix[j, ]))
        i = i + 1
    }
    
    Time = X[, 1]
    par(mfrow = c(1, 2))
    
    #png(filename = "timeVsPopulation", width = 800, height = 600) save plots
    
    plot(Time, X[, 2], ylab="S", ylim=c(0, max(X[nrow(X), 2:3])), type="l", lwd=2, main=paste("Susceptible Population vs Time for Starting Size", sampleSize), cex.main=0.70)
    plot(Time, X[, 3], ylab="R", ylim=c(0, max(X[nrow(X), 2:3])), type="l", lwd=2, main=paste("Resistant Population vs Time for Starting Size", sampleSize), cex.main=0.70)
    
    write.csv(X, file=paste("/Users/andrewhsu/Projects/PREP-NURA/NURA/data/cellData", sampleSize, ".birth_", l, ".death_", d, ".transformation_", h, ".mutation_", u1, ".csv", sep=""), row.names=FALSE)
}

