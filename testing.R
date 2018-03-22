library(slimr)
vcf <- read_vcf("~/projects/nea-over-time/data/simulations/deserts_tf_binding_site_h_0.5_rep_2_gen_2200.vcf.gz")

mut_info(vcf)
mut_info(vcf, mut_type = 0)
mut_info(vcf, mut_type = 1)
mut_info(vcf, pop_origin = 1)
mut_info(vcf, pop_origin = 2)
mut_info(vcf, t_min = 60000, t_max = 61000)

