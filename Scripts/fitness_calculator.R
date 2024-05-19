# Data preparation ====
library(rio)
library(tidyverse)

gene_table <- import("Carbon_source_wmy//4518-new_gene_table.xlsx", sheet = 1)
for_genes <- gene_table %>% 
  select(scaffoldId, begin, locusId)

strain_info <- import("Tn-Seq-blastn-short/all_blastn_short/all_blastn_unique_in_central_barcode.xlsx")

countT0 <- read_delim("Carbon_source_wmy/result/nb_count_bc.txt")
countCond <- read_delim("Carbon_source_wmy/result/ara_count_bc.txt")

count_merge <- merge(countCond, countT0, by = "barcode", all = TRUE)
colnames(count_merge)[2] <- "countCond"
colnames(count_merge)[3] <- "countT0"
count_merge <- 
  count_merge %>% 
  mutate(countCond = replace_na(countCond, 0)) %>% 
  mutate(countT0 =  replace_na(countT0, 0))

get_count <- merge(strain_info, count_merge, by = "barcode")

get_count <- 
  get_count %>% 
  mutate(genesused12 = if_else(f < 0.5 , 1, 2))

colnames(get_count)[10] <- 'locusId'

GeneFitness(for_genes, get_count[c("locusId", "f")], get_count$countCond, get_count$countT0,
            strainsUsed = NULL, genesUsed = NULL, genesUsed12 = get_count$locusId)

# FEBA.R -- analysis scripts for barcode sequencing data=======
#
# Uses mclapply() from the R parallel package to analyze experiments in parallel --
# set MC_CORES to control the #CPUs used (default is 2).
library(parallel);

# The key routines are:
#
# FEBA_Fit() -- analyze many fitness experiments with AvgStrainFitness() and NormalizeByScaffold()
#      returns a complex data structure
# FEBA_Save_Tables() -- Save the fitness data structure to a mini web-site, tab-delimited files and an R image
#
# Also see:
# AvgStrainFitness -- compute fitness values from counts for the post-experiment counts
#     and the Time0 counts
# NormalizeByScaffold -- normalize the fitness values by scaffold and position
# GeneFitness -- combines the above two, and also computes a t-like test statistic ("t").
#	To do this, it also computes fitness values for the 1st and 2nd half of most genes
#
# Limitations:
# High memory usage (~10GB to process 200K strains x 500 experiments)


# GeneFitness():
# genes -- must include locusId, scaffoldId, and begin
# strainInfo -- must include locusId and (unless use1 is overridden) f, the fraction of the gene
# 	     that the insertion is at
# countCond and countT0 -- counts for each strain
# strainsUsed & genesUsed -- see AvgStrainFitness()
# genesUsed12 -- ditto, for 1st and 2nd half fitness values
# use1 -- which strains are in 1st half (regardless of whether they are usable or not)
# other arguments are passed on to AvgStrainFitness()
# base_se -- likely amount of error in excess of that given by variation within fitness values
# 	for strains in a gene, due to erorrs in normalization or bias in the estimator
#
# Returns a data frame with a row for each gene in genesUsed. It includes
# locusId,
# fit (unnormalized), fitnorm (normalized),
# fit1 or fit2 for 1st- or 2nd-half (unnormalized, may be NA),
# fitnorm1 or fitnorm2 for normalized versions,
# se (estimated standard error of measurement), and t (the test statistic),
# as well as some other values from AvgStrainFitness(), notably sdNaive,
# which is a different (best-case) estimate of the standard error.
GeneFitness = function(genes, strainInfo, countCond, countT0,
                       strainsUsed, genesUsed, genesUsed12,
                       use1 = strainInfo$f < 0.5,
                       base_se = 0.1,
                       minGenesPerScaffold=10,
                       ...) {
  d = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
                       strainsUsed=strainsUsed, genesUsed=genesUsed, ...);
  d$fitnorm = NormalizeByScaffold(d$fit, d$locusId, genes, minToUse=minGenesPerScaffold);
  
  d1 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
                        strainsUsed=strainsUsed & use1, genesUsed=genesUsed12, ...);
  d2 = AvgStrainFitness(countCond, countT0, strainInfo$locusId,
                        strainsUsed=strainsUsed & !use1, genesUsed=genesUsed12, ...);
  if(any(as.character(d1$locusId) != as.character(d2$locusId))) stop("Non-matching locusId");
  
  i = match(d$locusId, d1$locusId);
  
  d$fit1 = d1$fit[i];
  d$fit2 = d2$fit[i];
  d$fitnorm1 = d$fit1 + (d$fitnorm-d$fit);
  d$fitnorm2 = d$fit2 + (d$fitnorm-d$fit);
  d$tot1 = d1$tot[i];
  d$tot1_0 = d1$tot0[i];
  d$tot2 = d2$tot[i];
  d$tot2_0 = d2$tot0[i];
  
  # for low n, the estimated variance is driven by the overall variance, which can be estimated
  # from the median difference between 1st and 2nd halves via the assumptions
  # Var(fit) = Var((fit1+fit2)/2) ~= Var(fit1-fit2)/4
  # median abs(normal variable) = qnorm(0.75) * sigma = 0.67 * sigma
  # which leads to Var(fit) = Var(fit1-fit2)/4
  # = sigma12**2/4 = median abs diff**2 / (qnorm(0.75)*2)**2
  # The median difference is used because a few genes may have genuine biological differences
  # between the fitness of the two halves.
  # Furthermore, assuming that genes with more reads are less noisy, this
  # pseudovariance should be rescaled based on sdNaive**2
  #
  pseudovar_std = median(abs(d$fit1-d$fit2),na.rm=T)**2 / (2*qnorm(0.75))**2;
  d$pseudovar = pseudovar_std * (d$sdNaive / median(d$sdNaive[!is.na(d$fit1)]))**2;
  # given the variable weighting in sumsq, it is not intuitive that the degrees of freedom is still n-1
  # however, this is the result given the assumption that the weighting is the inverse of the variance
  est_var = (d$pseudovar + d$sumsq)/d$n;
  d$se = sqrt(est_var);
  d$t = d$fitnorm/sqrt(base_se**2 + pmax(d$sdNaive**2, est_var));
  return(d);
}

# AvgStrainFitness():
# strainCounts -- counts at the end of the experiment condition
# strainT0 -- counts for Time0 for each strain
# strainLocus -- which locus the strain is associated with, or NA
#
# If genesUsed (as a list of locusId) and strainsUsed (as boolean vector) are provided,
# then considers only those strains & genes; minimum requirements.
#
# if returnStrainInfo is set, returns a list of two data frames, "genes" and "strains"
# normally returns a per-gene data frame.
#
# debug is for testing purposes.
AvgStrainFitness = function(strainCounts, strainT0, strainLocus,
                            minStrainT0 = 4, minGeneT0 = 40,
                            genesUsed=NULL, strainsUsed=NULL,
                            # maxWeight of N corresponds to having N reads on each side (if perfectly balanced); use 0 for even weighting
                            # 20 on each side corresponds to a standard error of ~0.5; keep maxWeight low because outlier strains
                            # often have higher weights otherwise.
                            # Use maxWeight = 0 for unweighted averages.
                            maxWeight = 20,
                            minGeneFactorNStrains=3,
                            returnStrainInfo=FALSE,
                            debug=getenv_numeric_or_default("FEBA_DEBUG",FALSE)) {
  if (length(strainCounts) < 1 || length(strainT0) < 1 || length(strainLocus) < 1
      || length(strainCounts) != length(strainT0) || length(strainCounts) != length(strainLocus))
    stop("No or misaligned input data");
  
  if (is.null(strainsUsed)) strainsUsed = strainT0 >= minStrainT0;
  if (is.null(genesUsed)) {
    geneT0 = aggregate(strainT0[strainsUsed], list(locusId=strainLocus[strainsUsed]), sum);
    genesUsed = geneT0$locusId[geneT0$x >= minGeneT0];
  }
  strainsUsed = strainsUsed & strainLocus %in% genesUsed;
  if (!any(strainsUsed)) stop("No usable strains");
  
  strainT0 = strainT0[strainsUsed];
  strainCounts = strainCounts[strainsUsed];
  strainLocus = strainLocus[strainsUsed];
  
  readratio = sum(strainCounts) / sum(strainT0);
  # use sqrt(readratio), or its inverse, instead of 1, so that the expectation
  # is about the same regardless of how well sampled the strain or gene is
  strainFit = mednorm(log2(sqrt(readratio) + strainCounts) - log2(1/sqrt(readratio) + strainT0));
  strainFitAdjust = 0;
  
  # Per-strain "smart" pseudocount to give a less biased per-strain fitness estimate.
  # This is the expected reads ratio, given data for the gene as a whole
  # Arguably, this should be weighted by T0 reads, but right now it isn't.
  # Also, do not do if we have just 1 or 2 strains, as it would just amplify noise
  
  # note use of as.vector() to remove names -- necessary for speed
  strainLocusF = as.factor(strainLocus);
  nStrains = table(strainLocusF);
  if(!all(names(nStrains)==levels(strainLocusF))) stop("Strain mismatch");
  nStrains = as.vector(nStrains);
  geneFit1 = mednorm(as.vector(tapply(strainFit, strainLocusF, median))); # used to use mean
  i = as.integer(strainLocusF); # from strain index to gene index
  strainPseudoCount = ifelse(nStrains[i] >= minGeneFactorNStrains, 2**geneFit1[i] * readratio, readratio);
  
  # And apportion the pseudocount equally (in log space) between condition-count and strain-count
  # to minimize the deviations from pseudocount = 1
  condPseudoCount = sqrt(strainPseudoCount);
  t0PseudoCount = 1/sqrt(strainPseudoCount);
  # (or could do some sort of weighted likelihood-based inference of fitness values, might be better)
  
  # for each strain: fitness, s.d., and weight
  strainFit = log2(condPseudoCount + strainCounts) - log2(t0PseudoCount + strainT0) - strainFitAdjust;
  strainSd = sqrt(1/(1+strainT0) + 1/(1+strainCounts)) / log(2);
  # use harmonic mean for weighting; add as small number to allow maxWeight = 0.
  strainWeight = 0.5 + pmin(maxWeight, 2/( 1/(1+strainT0) + 1/(1+strainCounts) ) );
  
  fitness = lapply(split(1:length(strainT0), list(locusId=strainLocus)),
                   function(j) {
                     n = length(j);
                     totw = sum(strainWeight[j]);
                     fitRaw = sum(strainWeight[j] * strainFit[j]) / totw;
                     tot = sum(strainCounts[j]);
                     tot0 = sum(strainT0[j]);
                     sd = sqrt(sum(strainWeight[j]**2 * strainSd[j]))/totw;
                     sumsq = sum(strainWeight[j] * (strainFit[j]-fitRaw)**2)/totw;
                     # high-N estimate of the noise in the log2 ratio of fitNaive
                     # But sdNaive is actually pretty accurate for small n -- e.g.
                     # simulations with E=10 on each side gave slightly light tails
                     # (r.m.s.(z) = 0.94).
                     sdNaive = sqrt( 1/(1+tot) + 1/(1+tot0) ) / log(2);
                     nEff = totw/max(strainWeight[j]);
                     c(fitRaw=fitRaw, sd=sd, sumsq=sumsq, sdNaive=sdNaive, n=n, nEff=nEff,
                       tot=tot, tot0=tot0);
                   });
  fitness = data.frame(do.call(rbind, fitness));
  fitness$fit = mednorm(fitness$fit);
  fitness$fitNaive = mednorm(log2(1+fitness$tot) - log2(1+fitness$tot0));
  fitness$locusId = row.names(fitness);
  if (is.integer(strainLocus)) fitness$locusId = as.integer(as.character(fitness$locusId));
  
  if(returnStrainInfo) return(list(genes=fitness,
                                   strains=data.frame(strainLocusF,strainCounts,strainT0,strainPseudoCount,strainFit,strainSd,strainWeight)));
  # else
  return(fitness);
}

# NormalizeByScaffold():
# values -- fitness values (as a vector)
# locusId -- the corresponding locusIds
# genes contains locusId, scaffoldId, and begin
# window -- window size for smoothing by medians. Must be odd, default 251. For scaffolds
#     with fewer genes than this, just uses the median.
# minToUse -- if a scaffold has too few genes, cannot correct for possible DNA extraction bias
# 	   so need to remove data for that gene (i.e., returns NA for them).
# returns
# normalized -- data with scaffold and position effects removed
#
NormalizeByScaffold = function(values, locusId, genes, window=251, minToUse=10, debug=FALSE) {
  i = match(locusId, genes$locusId);
  if(any(is.na(i))) stop("Fitness data for loci not in genes");
  beg = genes$begin[i];
  
  perScaffoldRows = split(1:length(values), genes$scaffoldId[i]);
  for (scaffoldId in names(perScaffoldRows)) {
    rows = perScaffoldRows[[scaffoldId]];
    if (length(rows) < minToUse) {
      if(debug) cat("Removing ",length(rows)," values for ", scaffoldId, "\n");
      values[rows] = NA;
    } else {
      med = median(values[rows]);
      if(debug) cat("Subtract median for ", scaffoldId, " ", med, "\n");
      values[rows] = values[rows] - med;
      
      if (length(rows) >= window) {
        # Used to use lowess
        # d = lowess(beg[rows], values[rows]);
        # if(debug) cat("Subtract loess for ", scaffoldId, " max effect is ", diff(range(d$y)), "\n");
        # values[rows] = values[rows] - approx(d$x,d$y, xout=beg[rows], rule=2, ties="ordered")$y;
        
        o = order(beg[rows]);
        m = runmed(values[rows[o]], window, endrule="constant");
        if(debug) cat("Subtract smoothed median for ", scaffoldId, " max effect is ",diff(range(m)), "\n");
        values[rows[o]] = values[rows[o]] - m;
        
        d = density(values[rows]);
        mode = d$x[which.max(d$y)];
        if (debug) cat("Subtract mode for ", scaffoldId, " which is at ", mode, "\n");
        values[rows] = values[rows] - mode;
      }
    }
  }
  return(values);
}

# median-based normalization
mednorm = function(x) x - median(x)