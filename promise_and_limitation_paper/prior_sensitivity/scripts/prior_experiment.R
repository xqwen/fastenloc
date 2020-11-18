#Initialize prior information
pd = 1e-3
pg = 5e-5


#now to obtain wide range of plausible a1 values
a1 = seq(from = -6, to = 8, by = .05)
pph4_vec = numeric(length(a1))

#Generate SCP for different levels of evidence (q1 and q2)
generate = function(q1,q2,alpha1) {
  for (i in 1:length(alpha1)) {
    BF1 = (q1/(1-q1))*((1-pd)/pd)
    BF2 = (q2/(1-q2))*((1-pg)/pg)
    a0 = log(pg/(pd*exp(alpha1[i])+1-pd-pg))
    pd1g1 = pd*exp(a0+alpha1[i])/(1+exp(a0+alpha1[i]))
    pd0g0 = (1-pd)/(1+exp(a0))
    pd1g0 = pd/(1+exp(a0+alpha1[i]))
    pd0g1 = (1-pd)*exp(a0)/(1+exp(a0))
    pph4_vec[i] = (pd1g1*BF1*BF2) / (pd0g0 + pd1g0 * BF1 + pd0g1 * BF2 + pd1g1*BF1*BF2)
  }
  return(pph4_vec)
}
pdf("results/SCP_vs_a1_t.pdf")
plot(a1,pph4_vec,xlab = expression(alpha[1]), type = "n",ylab = "Posterior Probability of Colocalization", ylim = c(0,1))
lines(a1,generate(0.9,0.9,a1), col = "red")
lines(a1,generate(0.9,0.5,a1),col = "blue", lty = 2)
lines(a1,generate(0.5,0.9,a1), col = "green")
lines(a1,generate(0.5,0.5,a1),col = "purple")
lines(a1, generate(0.05,0.05,a1),col = "black")
legend("topleft",legend = c("Strong eQTL and Strong GWAS","Strong eQTL and Moderate GWAS","Moderate eQTL and Strong GWAS", "Moderate eQTL and Moderate GWAS", "Weak eQTL and Weak GWAS"),col = c("red","blue","green","purple","black"),lty = c(1,2,1,1,1),cex = 0.7)
dev.off()

