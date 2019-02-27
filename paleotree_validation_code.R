# R script used to compare the output of FossilSim and paleotree

library(FossilSim)
library(paleotree)

spec = 0.5
ext = 0.2
foss = 0.15
ana = 0.1
propb = 0.6
n_extant = 20

paleotree_sim_raw = list()
while (length(paleotree_sim_raw) < 10000) {
  tmp = try(simFossilRecord(spec, ext, foss, anag.rate = ana, prop.bifurc = propb, nruns = 10, nExtant = n_extant))
  if(class(tmp) != "try-error") paleotree_sim_raw = c(paleotree_sim_raw, tmp)
}

paleotree_sim = lapply(paleotree_sim_raw, function(x){
  FossilSim::paleotree.record.to.fossils(x)
})

fs_trees = TreeSim::sim.bd.taxa(n_extant, 10000, spec, ext, complete = T)
fossilsim_sim = lapply(fs_trees, function(tr) {
  tax = FossilSim::sim.taxonomy(tr, beta = propb, lambda.a = ana)
  fos = FossilSim::sim.fossils.poisson(foss, taxonomy = tax)
  list(tree = tr, fossils = fos, taxonomy = tax)
})

fs_stats = function(sim_ds) {
  list(
    root_time = sapply(sim_ds, function(f) {
      max(ape::node.depth.edgelength(f$tree)) + f$tree$root.edge
    }),
    nfossils = sapply(sim_ds, function(f) {
      length(f$fossils$sp)
    }),
    nspec = sapply(sim_ds, function(f) {
      length(unique(f$taxonomy$sp))
    })
  )
}

paleotree_stats = fs_stats(paleotree_sim)
fossilsim_stats = fs_stats(fossilsim_sim)

comp_df = data.frame(fossilsim_stats)
comp_df$sim = "FossilSim"
tmp = data.frame(paleotree_stats)
tmp$sim = "Paleotree"
comp_df = rbind(comp_df,tmp)

par(mfrow=c(1,3))
boxplot(root_time ~ sim, comp_df, main = "Root time of tree")
boxplot(nspec ~ sim, comp_df, main = "Number of species in taxonomy")
boxplot(nfossils ~ sim, comp_df, main = "Number of fossils sampled")
