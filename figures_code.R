
library(FossilSim)

# new figure 2

set.seed(123)

t = TreeSim::sim.bd.taxa(n = 4, numbsim = 1, lambda = 1, mu = 0.3)[[1]]

set.seed(111) # 111

f = sim.fossils.poisson(rate = 2, tree = t)

pdf("fig2a.pdf", width = 5, height = 4)
plot(f, tree = t)
dev.off()

pdf("fig2b.pdf", width = 5, height = 4)
plot(f, tree = t, show.fossils = FALSE, show.ranges = TRUE) # show.strata = TRUE, strata = 4
dev.off()

# new examples figure
set.seed(122)
t = TreeSim::sim.bd.taxa(n = 8, numbsim = 1, lambda = 1, mu = 0.3)[[1]]
s = sim.taxonomy(tree = t, beta = 0.5, lambda.a = 1)
f = sim.fossils.poisson(rate = 3, taxonomy = s)

pdf("fig3a.pdf", width = 6, height = 6)
plot(s, tree = t, legend.position = "bottomright")
dev.off()

pdf("fig3b.pdf", width = 6, height = 6)
plot(f, tree = t, taxonomy = s, show.taxonomy = TRUE, show.ranges = TRUE)
dev.off()

set.seed(112) # 112
t = TreeSim::sim.bd.taxa(n = 8, numbsim = 1, lambda = 1, mu = 0.3)[[1]]


set.seed(110) # 
rate = 1

# define the distribution used to sample new rates
# in this case an exponential with a mean ~ 3
dist = function() { rexp(1, 1/4) }

# define the probability of the trait value changing at each speciation event
change.pr = 0.5

# simulate trait values under the independent model
rates = sim.trait.values(rate, tree = t, model = "independent", dist = dist, change.pr = change.pr)

# simulate fossils
f = sim.fossils.poisson(t, rate = rates)

pdf("fig3c.pdf", width = 6, height = 6)
plot(f, t)
dev.off()

set.seed(128) # 125 128 110 103 108  / 117 129

t = FossilSim::sim.fbd.taxa(n = 10, numbsim = 1, lambda = 1, mu = 0.5, psi = 2, 
                   complete = TRUE)[[1]]

pdf("fig3d.pdf", width = 6, height = 6)
FossilSim::rangeplot.asymmetric(t, complete = TRUE)
dev.off()


