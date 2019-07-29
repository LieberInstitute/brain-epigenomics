tab = read_excel("./Dropbox/wgbs_development/Third submission/Figures/Supplementary Tables/Supplementary_tables.xlsx", sheet = 1, skip = 2, range = "A3:T78")
num = read.csv("./Desktop/UMR_LMR_number.csv")
num = num[num$rep=="Discovery" & num$category=="LMR" & num$celltype!="Prenatal",
          colnames(num)!="X"]
num$age = tab$Age[match(num$id, tab$`Data ID`)]
num$celltype = factor(num$celltype)

num = split(num, num$celltype)
num = lapply(num, function(x) lm(num ~ age, x))
lapply(num, coef)

#$Glia
#(Intercept)         age 
# 65421.8900    138.2346 

#$Neuron
#(Intercept)         age 
# 64874.7760   -617.8338 

(coef(num$Neuron)[2] - coef(num$Glia)[2]) / coef(num$Glia)[2]
# -5.469458


genome = read.csv("./Desktop/UMR_LMR_genomeCoverage.csv")
genome = genome[genome$rep=="Discovery" & genome$category=="LMR" &
                  genome$celltype!="Prenatal", colnames(genome)!="X"]
genome$age = tab$Age[match(genome$id, tab$`Data ID`)]
genome$celltype = factor(genome$celltype)

genome = split(genome, genome$celltype)
genome = lapply(genome, function(x) lm(percent ~ age, x))
lapply(genome, coef)

#$Glia
# (Intercept)         age 
# 1.319276842 0.001902723 

#$Neuron
#(Intercept)         age 
#  1.2097065  -0.0146946 

(coef(genome$Neuron)[2] - coef(genome$Glia)[2]) / coef(genome$Glia)[2]
# -8.722929

