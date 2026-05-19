counts <- matrix(round(runif(1000,1,100)), ncol=4)
colnames(counts) <- c("ctrl1","ctrl2","treat1","treat2")
write.table(counts, "count/counts.txt", sep="\t")
cat("定量完成\n")