brnd_anova <- function(resp, fact, s) {

	output <- NULL

	for (i in 1:length(resp)) {

		# Single Factor Randomization ANOVA.
		# Output is a vector with DF treatment, DF error, MS error, F statistic, and p value.
		Fstat <- summary(aov(resp[,i] ~ fact))[[1]]$F[1]
		MSe <- summary(aov(resp[,i] ~ fact))[[1]]$Mean[2]
		DFt <- summary(aov(resp[,i] ~ fact))[[1]]$Df[1]
		DFe <- summary(aov(resp[,i] ~ fact))[[1]]$Df[2]

		perm_resp <- rep(0,length(resp[,i])*s)
		dim(perm_resp) <- c(length(resp[,i]),s)
		for (j in 1:s) {
			perm_resp[,j] <- sample(resp[,i],length(resp[,i]),replace=FALSE)}

		Fnull <- rep(0,s)
		for (j in 1:s) {
			temp <- summary(aov(perm_resp[,j] ~ fact))
			Fnull[j] <- temp[[1]]$F[1]}

		pval <- sum(Fnull > Fstat) / length(Fnull)

		ranova <- c(DFt,DFe,MSe,Fstat,pval)
		names(ranova) <- c("DFt","DFe","MSe","F","p")

		output <- rbind(output,ranova)
		}

	rownames(output) <- colnames(resp)
	output
	}

#Example
#table1 <- read.table(file="duploss_histories_branchlengthnormalized")
#brnd_anova(table1[,2:19553], table1[,1], 1000) -> results
#write.table(results, file="duploss_histories_branchlengthnormalized_results")

