require(Fantasio)
# data
require(HGDP.CEPH)
filepath <- system.file("extdata", "hgdp_ceph.bed", package="HGDP.CEPH")
x <- read.bed.matrix(filepath)
x <- set.stats(x)
x <- set.dist(x, HumanGeneticMap::genetic.map.b36)

x.be <- select.inds(x, population == "Bedouin")

# run by Hotspots
set.seed(1)
F1 <- Fantasio(bedmatrix=x.be, segments="Hotspots", segment.options = list(hotspots = hotspot_hg18), n = 5, verbose = TRUE)

# run by Segments
F3 <- Fantasio(bedmatrix=x.be, segments="Distance", n = 5, verbose = TRUE)

