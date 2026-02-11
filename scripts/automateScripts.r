# optimizedModel.r
require(dplyr)
require(ggplot2)
require(patchwork)
source("/Users/andrewhsu/Projects/PREP-NURA/NURA/scripts/optimizedModel.r")

automateRuns = function(modelType, argToVary, valsToTest){
    # Variables
    dataDir = "/Users/andrewhsu/Projects/PREP-NURA/NURA/data/Automated-Simulations"
    plotDir = "/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/Automated-Simulations"
    
    combFunctionArgs = list(
        birthRate = 0.5,
        deathRate = 0.5,
        turnOverRatio = 0.5,
        mutationRate1 = 1e-3,
        mutationRate2 = 1e-3,
        mutationRateB = 1e-3,
        drugInducedDeathRate1 = 0.01,
        drugInducedDeathRate2 = 0.01,
        sSize = 1000,
        r1Size = 0,
        r2Size = 0,
        rbSize = 0,
        maxTime = 100,
        tau = 0.1
    )
    
    cycFunctionArgs = list(
        birthRate = 0.5,
        deathRate = 0.5,
        turnOverRatio = 0.5,
        mutationRate1 = 1e-3,
        mutationRate2 = 1e-3,
        mutationRateB = 1e-3,
        drugInducedDeathRate1 = 0.01,
        drugInducedDeathRate2 = 0.01,
        sSize = 1000,
        r1Size = 0,
        r2Size = 0,
        rbSize = 0,
        tau = 0.1,
        cycleTime1 = 5,
        cycleTime2 = 10,
        totalCycles = 100
    )
    
    # File Naming
    fileName = NULL
    if(modelType == "Combinational"){
        fileName = paste("birthrate.", combFunctionArgs$birthrate, "_deathrate.", combFunctionArgs$deathrate, "_mutationrate1.", combFunctionArgs$mutationrate1, "_mutationrate2.", combFunctionArgs$mutationrate2, "_mutationrateb.", combFunctionArgs$mutationrateb, "_druginduceddeathrate1.", combFunctionArgs$druginduceddeathrate1, "_druginduceddeathrate2.", combFunctionArgs$druginduceddeathrate2, "_sSize.", combFunctionArgs$sSize, "_r1Size.", combFunctionArgs$r1Size, "_r2Size.", combFunctionArgs$r2Size, "_rbSize.", combFunctionArgs$rbSize, "_tau.", combFunctionArgs$tau, sep="")
    }
    else{
        fileName = paste("birthrate.", cycFunctionArgs$birthrate, "_deathrate.", cycFunctionArgs$deathrate, "_mutationrate1.", cycFunctionArgs$mutationrate1, "_mutationrate2.", cycFunctionArgs$mutationrate2, "_mutationrateb.", cycFunctionArgs$mutationrateb, "_druginduceddeathrate1.", cycFunctionArgs$druginduceddeathrate1, "_druginduceddeathrate2.", cycFunctionArgs$druginduceddeathrate2, "_sSize.", cycFunctionArgs$samplesize, "_r1Size.", cycFunctionArgs$r1Size, "_r2Size.", cycFunctionArgs$r2Size, "_rbSize.", cycFunctionArgs$rbSize, "_tau.", cycFunctionArgs$tau, sep="")
    }
    argToVaryIdx = regexpr(argToVary, fileName)[1]
    fileNameSubStr = substr(fileName, argToVaryIdx, nchar(fileName))
    underScoreIdx = regexpr("_", fileNameSubStr)[1] + argToVaryIdx - 1
    fileName = sub(substr(fileName, argToVaryIdx, underScoreIdx), "", fileName)
    
    if(argToVary == "turnOverRatio"){
        argToVaryIdx = regexpr("deathRate", fileName)[1]
        fileNameSubStr = substr(fileName, argToVaryIdx, nchar(fileName))
        underScoreIdx = regexpr("_", fileNameSubStr)[1] + argToVaryIdx - 1
        fileName = sub(substr(fileName, argToVaryIdx, underScoreIdx), "", fileName)
    }
    
    # Running
    set.seed(42)
    simulationRuns = lapply(valsToTest, function(x){
        currentArgs = NULL 
        
        if(argToVary == "turnOverRatio"){
            currentArgs$turnOverRatio = x
            currentArgs$deathRate = x * currentArgs$birthRate
        } 
        else{
            currentArgs[[argToVary]] = x
        }
        
        if(grepl("Combinational", modelType)){
            currentArgs = combFunctionArgs
            simulationResults = do.call(optimizedKomarova2D, currentArgs)
        } 
        else{
            currentArgs = cycFunctionArgs
            simulationResults = do.call(optimizedCyclicalTreatment, currentArgs)
        }
        
        summaryStats = as.data.frame(t(simulationResults[[2]]))
        summaryStats$ArgToVary = argToVary
        summaryStats$ArgValue = round(x, 5) 
        
        if(argToVary == "turnOverRatio"){
            simulationResults = simulationResults[[1]]
            simulationResults$turnOverRatio = round(currentArgs$turnOverRatio, 5)    
            simulationResults$deathRate = round(currentArgs$deathRate, 5) 
            simulationResults[[argToVary]] = round(x, 5)
        } 
        else{
            simulationResults = simulationResults[[1]]
            simulationResults[[argToVary]] = round(x, 5)
        }
        
        return(list(simulationResults, summaryStats))
    })
    simulationData = bind_rows(lapply(simulationRuns, "[[", 1))
    summaryData = bind_rows(lapply(simulationRuns, "[[", 2))
    names(summaryData) = c("Time", "TimeToSpawn", "MaxS", "MaxR1", "MaxR2", "MaxRB", "AvgS", "AvgR1", "AvgR2", "AvgRB", "VarS", "VarR1", "VarR2", "VarRB", "Success", "ArgToVary", "ArgValue")
    
    
    # Plotting
    subsetValsToTest = valsToTest[seq(1, length(valsToTest), length.out=5)]
    plotData = simulationData |> 
        filter(.data[[argToVary]] %in% subsetValsToTest) |>
        pivot_longer(cols=c("S", "R1", "R2", "RB"),
                     names_to="ResistanceType",
                     values_to="Count") |>
        mutate(ResistanceType=factor(ResistanceType, levels=c("S", "R1", "R2", "RB")))
    
    bigGridLog10 = ggplot(plotData, aes(x=Time, y=Count + 1, color=ResistanceType)) +
        geom_line(linewidth=1) +
        scale_y_log10() + 
        facet_grid(rows=vars(.data[[argToVary]]), cols=vars(ResistanceType), scales="free_y") +
        theme_bw() +
        labs(title=paste("Sensitivity Analysis (log10) Varying:", argToVary))
    
    bigGrid = ggplot(plotData, aes(x=Time, y=Count, color=ResistanceType)) +
        geom_line(linewidth=1) +
        facet_wrap(vars(.data[[argToVary]], ResistanceType), scales="free_y", ncol=4) +
        theme_bw() +
        labs(title=paste("Sensitivity Analysis Varying:", argToVary))
    
    
    # Saving
    saveDataPath = paste(dataDir, modelType, argToVary, sep="/")
    savePlotPath = paste(plotDir, modelType, argToVary, sep="/")
    
    dir.create(saveDataPath, recursive=TRUE, showWarnings=FALSE)
    dir.create(savePlotPath, recursive=TRUE, showWarnings=FALSE)
    
    write.csv(simulationData, file=paste(saveDataPath, "/", fileName, ".data.csv", sep=""), row.names=FALSE)
    write.csv(summaryData, file=paste(saveDataPath, "/", fileName, ".summary.csv", sep=""), row.names=FALSE)
    ggsave(paste(savePlotPath, "/", fileName, ".plotLog10.png", sep=""), bigGridLog10, width=14, height=10, units="in")
}

main = function(){
    modelType = c("Combinational", "Cyclical")
    combArgList = list(
        deathRate = c(seq(0.5, 0.9, 0.05), seq(0.9, 1, 0.01), seq(1, 2, 0.1)),
        turnOverRatio = seq(0.5, 1, 0.1),
        mutationRate1 = 10^seq(-9, -1, 2),
        mutationRate2 = 10^seq(-9, -1, 2),
        mutationRateB = 10^seq(-9, -1, 2),
        drugInducedDeathRate1 = 10^seq(-5, -1, 1),
        drugInducedDeathRate2 = 10^seq(-5, -1, 1),
        sSize = c(1000, 250, 500, 500, 0, 0, 0),
        r1Size = c(0, 250, 500, 0, 500, 1000, 0),
        r2Size = c(0, 250, 0, 500, 500, 0, 1000),
        rbSize = c(0, 250, 0, 0, 0, 0, 0),
        tau = 10^seq(-4, -1, 1)
    )
    
    cycArgList = list(
        deathRate = c(seq(0.5, 0.9, 0.05), seq(0.9, 1, 0.01), seq(1, 2, 0.1)),
        turnOverRatio = seq(0.5, 1, 0.1),
        mutationRate1 = 10^seq(-9, -1, 2),
        mutationRate2 = 10^seq(-9, -1, 2),
        mutationRateB = 10^seq(-9, -1, 2),
        drugInducedDeathRate1 = 10^seq(-5, -1, 1),
        drugInducedDeathRate2 = 10^seq(-5, -1, 1),
        sSize = c(1000, 250, 500, 500, 0, 0, 0),
        r1Size = c(0, 250, 500, 0, 500, 1000, 0),
        r2Size = c(0, 250, 0, 500, 500, 0, 1000),
        rbSize = c(0, 250, 0, 0, 0, 0, 0),
        tau = 10^seq(-4, -1, 1),
        cycleTime1 = c(5, 10),
        cycleTime2 = c(10, 5)
    )
    
    for(model in modelType){
        if(model == "Combinational"){
            argList = combArgList
        }
        else{
            argList = cycArgList
        }
        
        for(arg in names(argList)){
            print(paste("Running", model, "and Varying", arg))
            automateRuns(model, arg, argList[[arg]])
        }
    }
    print("Simulations complete.")
}

main()
