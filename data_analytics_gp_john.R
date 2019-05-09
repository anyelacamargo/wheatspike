
source('genpred_library.R')
#get_data()

load('w8data')
names(d)
# trait_data <- d$traitdata
# file <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip"
# DOex <- read_cross2(file)

cl <- makeCluster(detectCores() -1)
registerDoParallel(cl)

magicex <- read_cross2('MAGIC/magic1.json')

pr <- calc_genoprob(magicex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)
k <- calc_kinship(apr, "loco")
out <- scan1(apr, magicex$pheno[,5])
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out, magicex$gmap)
