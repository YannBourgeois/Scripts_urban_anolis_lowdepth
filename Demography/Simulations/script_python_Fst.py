import allel
import numpy as np

g = allel.GenotypeArray(allel.read_vcf("simulation_neutral.vcf",fields='calldata/GT')['calldata/GT'])
subpops = [range(0,32),range(32,44)]
ac1 = g.count_alleles(subpop=subpops[0])
ac2 = g.count_alleles(subpop=subpops[1])
num, den = allel.hudson_fst(ac1, ac2)
num

den

####Fst for the window
fst = np.sum(num) / np.sum(den)
singletons=allel.joint_sfs(ac1[:, 1],ac2[:, 1])[0,1] + allel.joint_sfs(ac1[:, 1],ac2[:, 1])[1,0]
print(fst,singletons,len(ac2))

