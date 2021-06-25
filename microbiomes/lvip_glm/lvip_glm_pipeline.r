## get user-edited environmental variables output_prefix, phy_fp, counts_fp, and model_matrix_fp
model_dir <- getwd()
source(file.path(model_dir, 'lvip_glm.env'))
##

## set up output directory
dir.create(output_prefix)
outfile <- file(file.path(output_prefix, 'runlog.log'), open = 'wt')
sink(outfile, type = 'output', split = TRUE)
##

counts <- as.matrix(read.table(counts_fp, sep='\t'))
counts <- counts[,colSums(counts) > 0]
modelMat <- read.table(model_matrix_fp, sep='\t')
factorNames <- as.character(modelMat[1,])
modelMat <- modelMat[-1,]
idx <- sapply(factorNames, function(x) which(x == unique(factorNames)))

library(ape)
microbeTree <- phytools::midpoint.root(read.tree(phy_fp))
finalMicrobeTree <- reorder(drop.tip(microbeTree, microbeTree$tip.label[!microbeTree$tip.label %in% colnames(counts)]), 'pruningwise')
finalMicrobeTree$edge.length <- finalMicrobeTree$edge.length / mean(finalMicrobeTree$edge.length)

counts <- counts[,finalMicrobeTree$tip.label]

## generate some summary numbers regarding microbes
microbeTips <- colnames(counts)
NT <- length(microbeTips)
NI <- finalMicrobeTree$Nnode
NN <- NI + NT
sa <- rbind(finalMicrobeTree$edge,c(0,NT+1))[NN:1,]
divergence <- c(finalMicrobeTree$edge.length[1:length(finalMicrobeTree$tip.label)],1,finalMicrobeTree$edge.length[(length(finalMicrobeTree$tip.label)+1):length(finalMicrobeTree$edge.length)])

##
NS <- nrow(modelMat)
NSB <- length(unique(factorNames))
NB_s <- ncol(modelMat)
##

relabund <- diag(1/rowSums(counts)) %*% counts
inv_log_max_contam <- 1 / log(max(apply(relabund, 2, function(x) max(x[x>0]) / min(x[x>0]))))


standat <- list(NS                 = NS,
                NI                 = NI,
                NT                 = NT,
                NB_s               = NB_s,
                NSB                = NSB,
                idx                = idx,
                count              = t(counts),
                divergence         = divergence,
                self               = sa[,2],
                ancestor           = sa[,1],
                X_s                = modelMat,
                inv_log_max_contam = inv_log_max_contam,
                shape_gnorm        = 7)


save.image(file.path(output_prefix, 'setup.RData'))

sampling_command <- paste('./lvip_glm',
                          paste0('data file=', file.path(output_prefix, 'data.json')),
                          #paste0('init=', file.path(output_prefix, 'inits.json')),
                          'output',
                          paste0('file=', file.path(output_prefix, 'samples_sampling.csv')),
                          paste0('refresh=', 1),
                          'method=sample algorithm=hmc',
                          'stepsize=1',
                          'engine=nuts',
                          'max_depth=10',
                          'adapt t0=10',
                          'delta=0.8',
                          'kappa=0.75',
                          #'init_buffer=10',
                          #'window=2',
                          'term_buffer=50',
                          'num_warmup=1000',
                          'num_samples=1000',
                          ('opencl platform=0 device=0')[opencl],
                          sep=' ')

cmdstanr::write_stan_json(standat, file.path(output_prefix, 'data.json'))

setwd(cmdstanr::cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], 'STANCFLAGS="--include-paths=', model_dir, '" ', file.path(model_dir, 'lvip_glm')))

setwd(model_dir)
print(sampling_command)
print(date())
system(sampling_command)

stan.fit <- cmdstanr::read_cmdstan_csv(file.path(output_prefix, 'samples_sampling.csv'),
                                       format = 'draws_array')

save.image(file.path(output_prefix, 'results.RData'))

monteCarloP <- function(x, pn='p') {
  if(pn == 'n') {
    res <- (1 + sum(x >= 0, na.rm=TRUE)) / (1 + length(x))
  } else if(pn == 'p') {
    res <- (1 + sum(x <= 0, na.rm=TRUE)) / (1 + length(x))
  }
} # https://arxiv.org/pdf/1603.05766.pdf

stan.fit.draws <- stan.fit$post_warmup_draws

sigma_prevalence <- stan.fit.draws[,,grep('^sigma_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_prevalence <- apply(sigma_prevalence,3,median)

sigma_abundance <- stan.fit.draws[,,grep('^sigma_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_abundance <- apply(sigma_abundance,3,median)

sd_prevalence <- stan.fit.draws[,,grep('^sd_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sd_prevalence <- apply(sd_prevalence,3,median)

sd_abundance <- stan.fit.draws[,,grep('^sd_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sd_abundance <- apply(sd_abundance,3,median)

theta_prevalence <- stan.fit.draws[,,grep('^theta_prevalence',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
theta_prevalence <- apply(theta_prevalence,3,median)

theta_abundance <- stan.fit.draws[,,grep('^theta_abundance',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
theta_abundance <- apply(theta_abundance,3,median)

inv_log_less_contamination <- stan.fit.draws[,,grep('^inv_log_less_contamination',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
inv_log_less_contamination <- apply(inv_log_less_contamination,3,median)

contaminant_overdisp <- stan.fit.draws[,,grep('^contaminant_overdisp',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
contaminant_overdisp <- apply(contaminant_overdisp,3,median)


delta_prevalence <- stan.fit.draws[,,grep('^delta_prevalence\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]

pos1 <- apply(delta_prevalence,c(2,3),monteCarloP,pn='p')
possig <- pos1 < 0.05
neg1 <- apply(delta_prevalence,c(2,3),monteCarloP,pn='n')
negsig <- neg1 < 0.05
anysig <- possig | negsig
delta_prevalence_sig <- array(anysig,dim=c(NB_s,NN))

which(delta_prevalence_sig, arr.ind = T)


delta_prevalence_med <- apply(delta_prevalence,3,median)
delta_prevalence_med <- array(delta_prevalence,dim=c(NB_s,NN))


delta_prevalence <- array(delta_prevalence,dim=c(dim(delta_prevalence)[c(1,2)],NB_s,NN))
phy2 <- finalMicrobeTree
delta_prevalence_from_root <- array(dim=dim(delta_prevalence))
for(fact in 1:dim(delta_prevalence)[[3]]) {
  for(chain in 1:dim(delta_prevalence)[[2]]) {
    for(iter in 1:dim(delta_prevalence)[[1]]) {
      phy2$edge.length <- delta_prevalence[iter,chain,fact,-(NT+1)]
      temp <- phytools::nodeHeights(phy2)[,2]
      delta_prevalence_from_root[iter,chain,fact,] <- c(temp[1:NT],delta_prevalence[iter,chain,fact,NT+1],temp[(NT+1):(NN-1)])
    }
  }
}


pos1 <- apply(delta_prevalence_from_root,c(3,4),monteCarloP,pn='p')
possig <- pos1 < 0.1
neg1 <- apply(delta_prevalence_from_root,c(3,4),monteCarloP,pn='n')
negsig <- neg1 < 0.1
anysig <- possig | negsig
delta_prevalence_from_root_sig <- array(anysig,dim=c(NB_s,NN))

which(delta_prevalence_from_root_sig, arr.ind = T)




delta_abundance <- stan.fit.draws[,,grep('^delta_abundance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]

pos1 <- apply(delta_abundance,c(2,3),monteCarloP,pn='p')
possig <- pos1 < 0.05
neg1 <- apply(delta_abundance,c(2,3),monteCarloP,pn='n')
negsig <- neg1 < 0.05
anysig <- possig | negsig
delta_abundance_sig <- array(anysig,dim=c(NB_s,NN))

delta_abundance_med <- apply(delta_abundance,3,median)
delta_abundance_med <- array(delta_abundance_med,dim=c(NB_s,NN))

delta_abundance <- array(delta_abundance,dim=c(dim(delta_abundance)[c(1,2)],NB_s,NN))
phy2 <- finalMicrobeTree
delta_abundance_from_root <- array(dim=dim(delta_abundance))
for(fact in 1:dim(delta_abundance)[[3]]) {
  for(chain in 1:dim(delta_abundance)[[2]]) {
    for(iter in 1:dim(delta_abundance)[[1]]) {
      phy2$edge.length <- delta_abundance[iter,chain,fact,-(NT+1)]
      temp <- phytools::nodeHeights(phy2)[,2]
      delta_abundance_from_root[iter,chain,fact,] <- c(temp[1:NT],delta_abundance[iter,chain,fact,NT+1],temp[(NT+1):(NN-1)])
    }
  }
}


pos1 <- apply(delta_abundance_from_root,c(3,4),monteCarloP,pn='p')
possig <- pos1 < 0.1
neg1 <- apply(delta_abundance_from_root,c(3,4),monteCarloP,pn='n')
negsig <- neg1 < 0.1
anysig <- possig | negsig
delta_abundance_from_root_sig <- array(anysig,dim=c(NB_s,NN))

tax <- read.table(taxa_fp, sep='\t', row.names=1)
tax[finalMicrobeTree$tip.label[phangorn::Descendants(finalMicrobeTree,2835)[[1]]],]
