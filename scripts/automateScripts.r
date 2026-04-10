# Re-run with a different turnover ratio (0.9-1)
# optimizedModel.r
require(dplyr)
require(ggplot2)
require(patchwork)
require(tidyr)
source("/Users/andrewhsu/Projects/PREP-NURA/NURA/scripts/optimizedModel.r")

automateRuns = function(modelType, argToVary, valsToTest){
    # Variables
    dataDir = "/Users/andrewhsu/Projects/PREP-NURA/NURA/data/Automated-Simulations"
    plotDir = "/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/Automated-Simulations"
    
    combFunctionArgs = list(
        birthRate = 0.5,
        deathRate = 0.5,
        turnOverRatio = 0.5,
        mutationRate1 = 1e-5,
        mutationRate2 = 1e-5,
        mutationRateB = 1e-5,
        drugInducedDeathRate1 = 0.001,
        drugInducedDeathRate2 = 0.001,
        sSize = 100,
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
        mutationRate1 = 1e-5,
        mutationRate2 = 1e-5,
        mutationRateB = 1e-5,
        drugInducedDeathRate1 = 0.001,
        drugInducedDeathRate2 = 0.001,
        sSize = 100,
        r1Size = 0,
        r2Size = 0,
        rbSize = 0,
        tau = 0.1,
        cycleTime1 = 5,
        cycleTime2 = 10,
        totalCycles = 10
    )
    
    # File Naming
    fileName = NULL
    if(modelType == "Combinational"){
        fileName = paste("birthrate.", combFunctionArgs$birthRate, "_deathrate.", combFunctionArgs$deathRate, "_mutationrate1.", combFunctionArgs$mutationRate1, "_mutationrate2.", combFunctionArgs$mutationRate2, "_mutationrateb.", combFunctionArgs$mutationRateB, "_druginduceddeathrate1.", combFunctionArgs$drugInducedDeathRate1, "_druginduceddeathrate2.", combFunctionArgs$drugInducedDeathRate2, "_sSize.", combFunctionArgs$sSize, "_r1Size.", combFunctionArgs$r1Size, "_r2Size.", combFunctionArgs$r2Size, "_rbSize.", combFunctionArgs$rbSize, "_tau.", combFunctionArgs$tau, sep="")
    }
    else{
        fileName = paste("birthrate.", cycFunctionArgs$birthRate, "_deathrate.", cycFunctionArgs$deathRate, "_mutationrate1.", cycFunctionArgs$mutationRate1, "_mutationrate2.", cycFunctionArgs$mutationRate2, "_mutationrateb.", cycFunctionArgs$mutationRateB, "_druginduceddeathrate1.", cycFunctionArgs$drugInducedDeathRate1, "_druginduceddeathrate2.", cycFunctionArgs$drugInducedDeathRate2, "_sSize.", cycFunctionArgs$sSize, "_r1Size.", cycFunctionArgs$r1Size, "_r2Size.", cycFunctionArgs$r2Size, "_rbSize.", cycFunctionArgs$rbSize, "_tau.", cycFunctionArgs$tau, sep="")
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
        if(grepl("Combinational", modelType)){
            currentArgs = combFunctionArgs
            simulationResults = do.call(optimizedKomarova2D, currentArgs)
        } 
        else{
            currentArgs = cycFunctionArgs
            simulationResults = do.call(optimizedCyclicalTreatment, currentArgs)
        }
        
        if(argToVary == "turnOverRatio"){
            currentArgs$turnOverRatio = x
            currentArgs$deathRate = x * currentArgs$birthRate
        } 
        else{
            currentArgs[[argToVary]] = x
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
    names(summaryData) = c("Time", "TimeToSpawn", "MaxS", "MaxR1", "MaxR2", "MaxRB", "AvgS", "AvgR1", "AvgR2", "AvgRB", "VarS", "VarR1", "VarR2", "VarRB", "Success", "ArgToVary", "ArgValue", "ModelType")
    
    
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
    saveDataPath = file.path(dataDir, modelType, argToVary)
    savePlotPath = file.path(plotDir, modelType, argToVary)
    
    dataFile = paste0(fileName, "_data.csv")
    summaryFile = paste0(fileName, "_summary.csv")
    plotFile = paste0(fileName, "_plotLog10.png")
    
    dir.create(saveDataPath, recursive=TRUE, showWarnings=FALSE)
    dir.create(savePlotPath, recursive=TRUE, showWarnings=FALSE)
    
    write.csv(simulationData, file=file.path(saveDataPath, dataFile), row.names=FALSE)
    write.csv(summaryData, file=file.path(saveDataPath, summaryFile), row.names=FALSE)
    ggsave(filename=file.path(savePlotPath, plotFile), bigGridLog10, width=14, height=10, units="in")
}

evaluateResistanceFactors = function(){
    print("Running Experiment: Determinants of Resistance")
    
    tumorSizes = c(10, 100, 1000)
    mutationRates = c(1e-5, 1e-4, 1e-3)
    drugTherapy = c("One Drug", "Two Drugs")
    results = data.frame()
    
    for(size in tumorSizes){
        for(rate in mutationRates){
            for(therapy in drugTherapy){
                for(rep in 1:10){
                    d1 = 0.01
                    d2 = ifelse(therapy == "Two Drugs", 0.01, 0.0)
                    sim = optimizedKomarova2D(sSize=size, mutationRate1=rate, mutationRate2=rate, mutationRateB=rate, drugInducedDeathRate1=d1, drugInducedDeathRate2=d2)
                    simData = sim[[1]]
                    rbEmerged = any(simData$RB > 0)
                    timeToRB = ifelse(rbEmerged, min(simData$Time[simData$RB > 0]), NA)
                    finalTumorSize = sum(tail(simData[, c("S", "R1", "R2", "RB")], 1))
                    
                    results = bind_rows(results, data.frame(
                        TumorSize = as.factor(size),
                        MutationRate = as.factor(rate),
                        Therapy = therapy,
                        RBEmerged = rbEmerged,
                        TimeToRB = timeToRB,
                        TreatmentFailed = finalTumorSize > size
                    ))
                }
            }
        }
    }
    
    summaryPlot = results |>
        group_by(TumorSize, MutationRate, Therapy) |>
        summarise(FailureRate = mean(TreatmentFailed), .groups="drop") |>
        ggplot(aes(x=TumorSize, y=FailureRate, fill=Therapy)) +
        geom_bar(stat="identity", position="dodge") +
        facet_wrap(~MutationRate, labeller=label_both) +
        theme_bw() +
        labs(title="Probability of Treatment Failure",
             y="Failure Probability",
             x="Initial Tumor Size")
    ggsave("~/Projects/PREP-NURA/NURA/plots/resistanceDeterminants.png", summaryPlot, width=10, height=6)
    print("Experiment Complete.")
}

evaluateResistanceTiming = function(numPatients=50){
    print("Running Experiment: Pre-treatment vs Acquired Resistance")
    results = data.frame()
    
    for(patient in 1:numPatients){
        phase1 = optimizedKomarova2D(sSize=100, drugInducedDeathRate1=0, drugInducedDeathRate2=0, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, maxTime=30)
        p1Data = phase1[[1]]
        startTreatState = tail(p1Data, 1)
        
        totalCells = startTreatState$S + startTreatState$R1 + startTreatState$R2 + startTreatState$RB
        if(totalCells == 0){
            results = bind_rows(results, data.frame(
                Patient=patient,
                Verdict="3_No Resistance",
                FinalBurden=0,
                Success=TRUE
            ))
            next
        }
        
        preExistingRB = startTreatState$RB > 0
        phase2 = optimizedKomarova2D(sSize=startTreatState$S, r1Size=startTreatState$R1, r2Size=startTreatState$R2, rbSize=startTreatState$RB, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, maxTime=30)
        p2Data = phase2[[1]]
        endTreatState = tail(p2Data, 1)
        
        acquiredRB = (!preExistingRB) & (endTreatState$RB > 0)
        verdict = "3_No Resistance"
        if(preExistingRB) verdict = "1_Pre-existing (Before Treatment)"
        if(acquiredRB) verdict = "2_Acquired (During Treatment)"
        
        finalTumorBurden = endTreatState$S + endTreatState$R1 + endTreatState$R2 + endTreatState$RB
        success = finalTumorBurden < 100
        
        results = bind_rows(results, data.frame(
            Patient=patient,
            Verdict=verdict,
            FinalBurden=finalTumorBurden,
            Success=success
        ))
    }
    timePlot = ggplot(results, aes(x=Verdict, y=FinalBurden + 1, fill=Verdict)) +
        geom_boxplot(alpha=0.6) +
        geom_jitter(width=0.2, alpha=0.8) +
        scale_y_log10() +
        theme_bw() +
        theme(legend.position="none") +
        labs(title="Impact of Resistance Timing on Treatment Outcome",
             y="Final Tumor Cell Count (log10)",
             x="Verdict of RB Cells")
    ggsave("~/Projects/PREP-NURA/NURA/plots/resistanceTimingImpact30.png", timePlot, width=8, height=6)
    print("Experiment Complete.")
}

# Change treatment to start when first doubly-resistant cell starts are before it even starts.
evaluateResistanceTimingCyc = function(numPatients=50){
    rm(.Random.seed)
    print("Running Experiment: Pre-treatment vs Acquired Resistance")
    results = data.frame()
    
    phase1Threshold = 1
    for(patient in 1:numPatients){
        phase1 = optimizedCyclicalTreatment(sSize=100, drugInducedDeathRate1=0, drugInducedDeathRate2=0, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, totalCycle=6, cycleTime1=2, cycleTime2=2)
        p1Data = phase1[[1]]
        startTreatState = tail(p1Data, 1)
        
        totalCells = startTreatState$S + startTreatState$R1 + startTreatState$R2 + startTreatState$RB
        if(totalCells == 0){
            results = bind_rows(results, data.frame(
                Patient=patient,
                Verdict="3_No Resistance",
                FinalBurden=0,
                Success=TRUE
            ))
            next
        }
        
        preExistingRB = startTreatState$RB > 0
        phase2 = optimizedCyclicalTreatment(sSize=startTreatState$S, r1Size=startTreatState$R1, r2Size=startTreatState$R2, rbSize=startTreatState$RB, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, totalCycle=6, cycleTime1=2, cycleTime2=2)
        p2Data = phase2[[1]]
        endTreatState = tail(p2Data, 1)
        
        acquiredRB = (!preExistingRB) & (endTreatState$RB > 0)
        verdict = "3_No Resistance"
        if(preExistingRB) verdict = "1_Pre-existing (Before Treatment)"
        if(acquiredRB) verdict = "2_Acquired (During Treatment)"
        
        finalTumorBurden = endTreatState$S + endTreatState$R1 + endTreatState$R2 + endTreatState$RB
        success = finalTumorBurden < 100
        
        results = bind_rows(results, data.frame(
            Patient=patient,
            Verdict=verdict,
            FinalBurden=finalTumorBurden,
            Success=success
        ))
    }
    timePlot = ggplot(results, aes(x=Verdict, y=FinalBurden + 1, fill=Verdict)) +
        geom_boxplot(alpha=0.6) +
        geom_jitter(width=0.2, alpha=0.8) +
        scale_y_log10() +
        theme_bw() +
        theme(legend.position="none") +
        labs(title="Impact of Resistance Timing on Treatment Outcome",
             y="Final Tumor Cell Count (log10)",
             x="Verdict of RB Cells")
    ggsave("~/Projects/PREP-NURA/NURA/plots/resistanceTimingImpactTotalCyc6Cyc2.png", timePlot, width=8, height=6)
    print("Experiment Complete.")
    set.seed(42)
}

# Zoom in plot (inset) of turnover ratio from 0.8 to 1
turnOverBaoPlot = function(){
    print("Running Experiment: Turnover Ratio vs Number of Drugs")
    turnOverRatios = c(seq(0.1, 0.9, by=0.1), seq(0.8, 0.9, by=0.25))
    results = data.frame()
    
    for(tor in turnOverRatios){
        for(rep in 1:10){
            currentDeathRate = 1.0 * tor
            
            sim1 = optimizedKomarova2D(birthRate=1.0, deathRate=currentDeathRate, turnOverRatio=tor, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0, sSize=100, maxTime=25)
            maxRB1 = max(sim1[[1]]$RB)
            
            sim2 = optimizedKomarova2D(birthRate=1.0, deathRate=currentDeathRate, turnOverRatio=tor, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, sSize=100, maxTime=25)
            maxRB2 = max(sim2[[1]]$RB)
            
            results = bind_rows(results, data.frame(
                TurnOverRatio = tor,
                MaxRB = maxRB1,
                Therapy = "1 Drug"
            ))
            
            results = bind_rows(results, data.frame(
                TurnOverRatio = tor,
                MaxRB = maxRB2,
                Therapy = "2 Drugs"
            ))
        }
    }
    ToRPlot = ggplot(results, aes(x=TurnOverRatio, y=MaxRB + 1, color=Therapy)) +
        geom_point(alpha=0.4, position=position_jitter(width=0.02)) +
        geom_smooth(method="loess", se=FALSE, linewidth=1.2) +
        scale_y_log10() +
        theme_bw() +
        labs(title="Effect of Turnover Ratio on Resistance Emergence",
             x="Turnover Ratio",
             y="Maximum RB Population (log10)",
             color="Treatment Type")
    
    insetData = results |> filter(TurnOverRatio >= 0.8)
    insetPlot = ggplot(insetData, aes(x=TurnOverRatio, y=MaxRB + 1, color=Therapy)) +
        geom_point(alpha=0.6, position=position_jitter(width=0.005)) +
        geom_smooth(method="loess", se=FALSE, linewidth=1) +
        scale_y_log10() +
        theme_bw() +
        theme(legend.position="none",
              plot.title=element_text(size=10),
              axis.title=element_blank()) +
        labs(title="Zoom: 0.8 to 0.9")
    
    finalPlot = ToRPlot + inset_element(insetPlot, left=0.05, bottom=0.055, right=0.45, top=0.85)
    
    ggsave("~/Projects/PREP-NURA/NURA/plots/turnOverVSDrugs.png", finalPlot, width=8, height=6)
    print("Experiment Complete")
}

# Calculate variance
mutationBaoPlot = function(){
    print("Running Experiment: Mutation vs Number of Drugs")
    mutRates = 10^seq(-7, -3, by=1)
    results = data.frame()
    
    for(m in mutRates){
        for(rep in 1:10){
            sim1 = optimizedKomarova2D(mutationRate1=m, mutationRate2=m, mutationRateB=m, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0, sSize=100, maxTime=25)
            maxRB1 = max(sim1[[1]]$RB)
            
            sim2 = optimizedKomarova2D(mutationRate1=m, mutationRate2=m, mutationRateB=m, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, sSize=100, maxTime=25)
            maxRB2 = max(sim2[[1]]$RB)
            
            results = bind_rows(results, data.frame(
                MutationRate = m,
                MaxRB = maxRB1,
                Therapy = "1 Drug"
            ))
            
            results = bind_rows(results, data.frame(
                MutationRate = m,
                MaxRB = maxRB2,
                Therapy = "2 Drugs"
            ))
        }
    }
    
    summaryData = results |>
        group_by(MutationRate, Therapy) |>
        summarise(
            MeanRB = mean(MaxRB),
            SE_RB = sd(MaxRB) / sqrt(n()),
            .groups="drop"
        )
    
    mutPlot = ggplot() +
        geom_point(data=results, aes(x=MutationRate, y=MaxRB + 1, color=Therapy), alpha=0.2, position=position_jitter(width=0.05)) +
        geom_errorbar(data=summaryData, aes(x=MutationRate, ymin=pmax(1, (MeanRB - SE_RB) + 1), ymax=(MeanRB + SE_RB) + 1, color=Therapy), width=0.15, linewidth=0.8) +
        geom_line(data=summaryData, aes(x=MutationRate, y=MeanRB + 1, color=Therapy), size=1) +
        scale_x_log10() +
        scale_y_log10() +
        theme_bw() +
        labs(title="Effect of Mutation on Resistance Emergence",
             x="Mutation Rate (log10)",
             y="Maximum RB Population (log10)",
             color="Treatment Type")
    ggsave("~/Projects/PREP-NURA/NURA/plots/mutationVSDrugs.png", mutPlot, width=8, height=6)
    print("Experiment Complete")
}

# Increase baseline drug-induced death rate to try to knock down RB population
numberOfDrugsPlot = function(){
    print("Running Experiment: Number of Drugs")
    noDrugs = optimizedKomarova2D(drugInducedDeathRate1=0, drugInducedDeathRate2=0)
    oneDrug = optimizedKomarova2D(drugInducedDeathRate1=0.5, drugInducedDeathRate2=0)
    twoDrugsx1 = optimizedKomarova2D(drugInducedDeathRate1=0.5, drugInducedDeathRate2=0.5)
    twoDrugsx2 = optimizedKomarova2D(drugInducedDeathRate1=0.5, drugInducedDeathRate2=1)
    
    df_no = noDrugs[[1]][, c("Time", "RB")]
    df_no$Treatment = "0 Drugs (Control)"
    
    df_one = oneDrug[[1]][, c("Time", "RB")]
    df_one$Treatment = "1 Drug (d1=0.5)"
    
    df_two1 = twoDrugsx1[[1]][, c("Time", "RB")]
    df_two1$Treatment = "2 Drugs (d1=0.5, d2=0.5)"
    
    df_two2 = twoDrugsx2[[1]][, c("Time", "RB")]
    df_two2$Treatment = "2 Drugs (d1=0.5, d2=1)"
    
    plotData = bind_rows(df_no, df_one, df_two1, df_two2)
    
    plotData$Treatment = factor(plotData$Treatment, 
                                levels = c("0 Drugs (Control)", 
                                           "1 Drug (d1=0.5)", 
                                           "2 Drugs (d1=0.5, d2=0.5)", 
                                           "2 Drugs (d1=0.5, d2=1)"))
    
    drugsPlot = ggplot(plotData, aes(x=Time, y=RB + 1, color=Treatment)) +
        geom_line(linewidth=1.2) +
        scale_y_log10() + # Log scale is best for resistance populations
        theme_bw() +
        labs(title="Emergence of Double Resistance (RB) over Time",
             subtitle="Comparing Monotherapy vs Combinational Regimens",
             x="Time",
             y="RB Cell Population (log10)",
             color="Treatment Regimen")
    
    ggsave("~/Projects/PREP-NURA/NURA/plots/numberOfDrugsTimeSeriesDIDR0.05.png", drugsPlot, width=8, height=6)
    print("Number of Drugs Experiment Complete.")
}

numberOfDrugsPlotCyc = function(){
    print("Running Experiment: Number of Drugs")
    noDrugs = optimizedCyclicalTreatment(drugInducedDeathRate1=0, drugInducedDeathRate2=0, totalCycles=6, cycleTime1=2, cycleTime2=2)
    oneDrug = optimizedCyclicalTreatment(drugInducedDeathRate1=0.1, drugInducedDeathRate2=0, totalCycles=6, cycleTime1=2, cycleTime2=2)
    twoDrugsx1 = optimizedCyclicalTreatment(drugInducedDeathRate1=0.1, drugInducedDeathRate2=0.1, totalCycles=6, cycleTime1=2, cycleTime2=2)
    twoDrugsx2 = optimizedCyclicalTreatment(drugInducedDeathRate1=0.1, drugInducedDeathRate2=0.2, totalCycles=6, cycleTime1=2, cycleTime2=2)
    
    df_no = noDrugs[[1]][, c("Time", "RB")]
    df_no$Treatment = "0 Drugs (Control)"
    
    df_one = oneDrug[[1]][, c("Time", "RB")]
    df_one$Treatment = "1 Drug (d1=0.1)"
    
    df_two1 = twoDrugsx1[[1]][, c("Time", "RB")]
    df_two1$Treatment = "2 Drugs (d1=0.1, d2=0.1)"
    
    df_two2 = twoDrugsx2[[1]][, c("Time", "RB")]
    df_two2$Treatment = "2 Drugs (d1=0.1, d2=0.2)"
    
    plotData = bind_rows(df_no, df_one, df_two1, df_two2)
    
    plotData$Treatment = factor(plotData$Treatment, 
                                levels = c("0 Drugs (Control)", 
                                           "1 Drug (d1=0.1)", 
                                           "2 Drugs (d1=0.1, d2=0.1)", 
                                           "2 Drugs (d1=0.1, d2=0.2)"))
    
    drugsPlot = ggplot(plotData, aes(x=Time, y=RB + 1, color=Treatment)) +
        geom_line(linewidth=1.2) +
        scale_y_log10() + # Log scale is best for resistance populations
        theme_bw() +
        labs(title="Emergence of Double Resistance (RB) over Time",
             subtitle="Comparing Monotherapy vs Combinational Regimens",
             x="Time",
             y="RB Cell Population (log10)",
             color="Treatment Regimen")
    
    ggsave("~/Projects/PREP-NURA/NURA/plots/numberOfDrugsTimeSeriesDIDR0.1Cyc.png", drugsPlot, width=8, height=6)
    print("Number of Drugs Experiment Complete.")
}

turnOverBaoPlotCyc = function(){
    print("Running Experiment: Turnover Ratio vs Number of Drugs")
    turnOverRatios = c(seq(0.1, 0.9, by=0.1), seq(0.8, 0.9, by=0.25))
    results = data.frame()
    
    for(tor in turnOverRatios){
        for(rep in 1:10){
            currentDeathRate = 1.0 * tor
            
            sim1 = optimizedCyclicalTreatment(birthRate=1.0, deathRate=currentDeathRate, turnOverRatio=tor, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0, sSize=100, totalCycles=6, cycleTime1=2, cycleTime2=2)
            maxRB1 = max(sim1[[1]]$RB)
            
            sim2 = optimizedCyclicalTreatment(birthRate=1.0, deathRate=currentDeathRate, turnOverRatio=tor, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, sSize=100, totalCycles=6, cycleTime1=2, cycleTime2=2)
            maxRB2 = max(sim2[[1]]$RB)
            
            results = bind_rows(results, data.frame(
                TurnOverRatio = tor,
                MaxRB = maxRB1,
                Therapy = "1 Drug"
            ))
            
            results = bind_rows(results, data.frame(
                TurnOverRatio = tor,
                MaxRB = maxRB2,
                Therapy = "2 Drugs"
            ))
        }
    }
    ToRPlot = ggplot(results, aes(x=TurnOverRatio, y=MaxRB + 1, color=Therapy)) +
        geom_point(alpha=0.4, position=position_jitter(width=0.02)) +
        geom_smooth(method="loess", se=FALSE, linewidth=1.2) +
        scale_y_log10() +
        theme_bw() +
        labs(title="Effect of Turnover Ratio on Resistance Emergence",
             x="Turnover Ratio",
             y="Maximum RB Population (log10)",
             color="Treatment Type")
    
    ggsave("~/Projects/PREP-NURA/NURA/plots/turnOverVSDrugsCyc.png", ToRPlot, width=8, height=6)
    print("Experiment Complete")
}

# Calculate variance
mutationBaoPlotCyc = function(){
    print("Running Experiment: Mutation vs Number of Drugs")
    mutRates = 10^seq(-7, -3, by=1)
    results = data.frame()
    
    for(m in mutRates){
        for(rep in 1:10){
            sim1 = optimizedCyclicalTreatment(mutationRate1=m, mutationRate2=m, mutationRateB=m, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0, sSize=100, totalCycles=6, cycleTime1=2, cycleTime2=2)
            maxRB1 = max(sim1[[1]]$RB)
            
            sim2 = optimizedCyclicalTreatment(mutationRate1=m, mutationRate2=m, mutationRateB=m, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, sSize=100, totalCycles=6, cycleTime1=2, cycleTime2=2)
            maxRB2 = max(sim2[[1]]$RB)
            
            results = bind_rows(results, data.frame(
                MutationRate = m,
                MaxRB = maxRB1,
                Therapy = "1 Drug"
            ))
            
            results = bind_rows(results, data.frame(
                MutationRate = m,
                MaxRB = maxRB2,
                Therapy = "2 Drugs"
            ))
        }
    }
    
    summaryData = results |>
        group_by(MutationRate, Therapy) |>
        summarise(
            MeanRB = mean(MaxRB),
            SE_RB = sd(MaxRB) / sqrt(n()),
            .groups="drop"
        )
    
    mutPlot = ggplot() +
        geom_point(data=results, aes(x=MutationRate, y=MaxRB + 1, color=Therapy), alpha=0.2, position=position_jitter(width=0.05)) +
        geom_errorbar(data=summaryData, aes(x=MutationRate, ymin=pmax(1, (MeanRB - SE_RB) + 1), ymax=(MeanRB + SE_RB) + 1, color=Therapy), width=0.15, linewidth=0.8) +
        geom_line(data=summaryData, aes(x=MutationRate, y=MeanRB + 1, color=Therapy), size=1) +
        scale_x_log10() +
        scale_y_log10() +
        theme_bw() +
        labs(title="Effect of Mutation on Resistance Emergence",
             x="Mutation Rate (log10)",
             y="Maximum RB Population (log10)",
             color="Treatment Type")
    ggsave("~/Projects/PREP-NURA/NURA/plots/mutationVSDrugsCyc.png", mutPlot, width=8, height=6)
    print("Experiment Complete")
}

evaluateTumorHeterogeneity = function(){
    print("Running Experiment: Impact of Initial Tumor Heterogeneity")
    startingStates = list(
        list(name="All Susceptible (100S)", sSize=100, r1Size=0, r2Size=0, rbSize=0),
        list(name="Some Single Res (90S, 10R1)", sSize=90, r1Size=10, r2Size=0, rbSize=0),
        list(name="Mixed Single Res (90S, 5R1, 5R2)", sSize=90, r1Size=5, r2Size=5, rbSize=0),
        list(name="Some Double Res (95S, 5RB)", sSize=95, r1Size=0, r2Size=0, rbSize=5)
    )
    
    simulationData = data.frame()
    for(state in startingStates){
        sim = optimizedKomarova2D(sSize=state$sSize, r1Size=state$r1Size, r2Size=state$r2Size, rbSize=state$rbSize, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, maxTime=25)
        df = sim[[1]]
        df$InitialState = state$name
        simulationData = bind_rows(simulationData, df)
    }
    plotData = simulationData |>
        pivot_longer(cols=c("S", "R1", "R2", "RB"),
                     names_to="ResistanceType",
                     values_to="Count") |>
        mutate(ResistanceType=factor(ResistanceType, levels=c("S", "R1", "R2", "RB")))
    heterogeneityPlot = ggplot(plotData, aes(x=Time, y=Count + 1, color=ResistanceType)) +
        geom_line(linewidth=1.2) +
        scale_y_log10() +
        facet_grid(rows=vars(InitialState), cols=vars(ResistanceType), scale="free_y") +
        theme_bw() +
        theme(legend.position="none",
              strip.text.y=element_text(angle=0)) +
        labs(title="Impact of Tumor Heterogeneity on Clonal Evolution",
             y="Cell Count (log10)",
             x="Time")
    ggsave("~/Projects/PREP-NURA/NURA/plots/tumorHeterogeneity.png", heterogeneityPlot, width=8, height=6)
    print("Experiment Complete")
}

evaluateTumorHeterogeneityCyc = function(){
    print("Running Experiment: Impact of Initial Tumor Heterogeneity")
    startingStates = list(
        list(name="All Susceptible (100S)", sSize=100, r1Size=0, r2Size=0, rbSize=0),
        list(name="Some Single Res (90S, 10R1)", sSize=90, r1Size=10, r2Size=0, rbSize=0),
        list(name="Mixed Single Res (90S, 5R1, 5R2)", sSize=90, r1Size=5, r2Size=5, rbSize=0),
        list(name="Some Double Res (95S, 5RB)", sSize=95, r1Size=0, r2Size=0, rbSize=5)
    )
    
    simulationData = data.frame()
    for(state in startingStates){
        sim = optimizedCyclicalTreatment(sSize=state$sSize, r1Size=state$r1Size, r2Size=state$r2Size, rbSize=state$rbSize, drugInducedDeathRate1=0.01, drugInducedDeathRate2=0.01, mutationRate1=1e-5, mutationRate2=1e-5, mutationRateB=1e-5, totalCycles=6, cycleTime1=2, cycleTime2=2)
        df = sim[[1]]
        df$InitialState = state$name
        simulationData = bind_rows(simulationData, df)
    }
    plotData = simulationData |>
        pivot_longer(cols=c("S", "R1", "R2", "RB"),
                     names_to="ResistanceType",
                     values_to="Count") |>
        mutate(ResistanceType=factor(ResistanceType, levels=c("S", "R1", "R2", "RB")))
    heterogeneityPlot = ggplot(plotData, aes(x=Time, y=Count + 1, color=ResistanceType)) +
        geom_line(linewidth=1.2) +
        scale_y_log10() +
        facet_grid(rows=vars(InitialState), cols=vars(ResistanceType), scale="free_y") +
        theme_bw() +
        theme(legend.position="none",
              strip.text.y=element_text(angle=0)) +
        labs(title="Impact of Tumor Heterogeneity on Clonal Evolution",
             y="Cell Count (log10)",
             x="Time")
    ggsave("~/Projects/PREP-NURA/NURA/plots/tumorHeterogeneityCyc.png", heterogeneityPlot, width=8, height=6)
    print("Experiment Complete")
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
        sSize = c(100, 25, 50, 50, 0, 0, 0),
        r1Size = c(0, 25, 50, 0, 50, 100, 0),
        r2Size = c(0, 25, 0, 50, 50, 0, 100),
        rbSize = c(0, 25, 0, 0, 0, 0, 0),
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
        sSize = c(100, 25, 50, 50, 0, 0, 0),
        r1Size = c(0, 25, 50, 0, 50, 100, 0),
        r2Size = c(0, 25, 0, 50, 50, 0, 100),
        rbSize = c(0, 25, 0, 0, 0, 0, 0),
        tau = 10^seq(-4, -1, 1),
        cycleTime1 = c(5, 10),
        cycleTime2 = c(10, 5)
    )
    
    # for(model in modelType){
    #     if(model == "Combinational"){
    #         argList = combArgList
    #     }
    #     else{
    #         argList = cycArgList
    #     }
    # 
    #     for(arg in names(argList)){
    #         print(paste("Running", model, "and Varying", arg))
    #         automateRuns(model, arg, argList[[arg]])
    #     }
    # }
    # evaluateResistanceFactors()
    # evaluateResistanceTiming(100)
    # evaluateResistanceTimingCyc(100)
    # turnOverBaoPlot()
    # turnOverBaoPlotCyc()
    # mutationBaoPlot()
    # mutationBaoPlotCyc()
    # numberOfDrugsPlot()
    # numberOfDrugsPlotCyc()
    # evaluateTumorHeterogeneity()
    evaluateTumorHeterogeneityCyc()
    
    print("Simulations complete.")
}

main()
