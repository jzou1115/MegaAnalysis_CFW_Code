library(ggplot2)
library(mvtnorm)
library(pracma)

args = commandArgs(trailingOnly=TRUE)
filename<-args[1]
threshold <- as.numeric(args[3])
correction <- args[4]

data<-read.table(filename, header=T)




#estimating variance of s1, accounting for missing data
MLE<-function(var, s1, sampleSizeS1, M, MSig, ZScore){
	estimate<-0
	for(i in s1){
		estimate<-estimate+log(dnorm(i, mean=0, sd=sqrt(var)))
	}
	estimate<-estimate+(M-MSig)*log(pnorm(ZScore*sqrt(sampleSizeS1), mean=0, sd=sqrt(var))-pnorm((-ZScore)*sqrt(sampleSizeS1), mean=0, sd=sqrt(var)))

	return(estimate)
}

#variance explained by genetics
expected_mean_ratio<-function(sigma_g_2, maxVar, n1){
	sigma_g_2/maxVar
}

#estimate variance components
estimate_parameters <- function(s1, s2, n1, n2, M, MSig, Zscore){
	#sigma_g, sigma_c1, sigma_c2 initial values
	parameters <- c(0.0000000000000000000000000001, 0.0000000000000000000000000001, 0.0000000000000000000000000001)


	######### estimate total variance in s1  #############
	max<- -Inf
	maxVar<- -Inf
	for(i in seq(from=0.1, to=10, by=.01)){
		temp<-MLE(i, s1, n1, M, MSig, Zscore)
		if (temp>max) {
			max<-temp
 			maxVar<-i
 		}
	}


	print(paste0("maxVar: ", toString(maxVar)))

	######### estimate sigma_g  #############
#	sigma_g_estimator <- 0.0000000000000000000000000001
	min_rms <- Inf
	var_g_est2 <- 0.0000000000000000000000000001
	#for (i in 1:100000){
	for(i in seq(from=0.0000000000000000000000000001,to=10, by=0.01)){
		ratio<-expected_mean_ratio(i, maxVar, n1)
		expected_s2 <- s1*ratio
  		cur_rms <- sqrt(sum((expected_s2-s2)^2))
  		if (cur_rms < min_rms){
    			min_rms<- cur_rms
    			var_g_est2 <- i
  		}
	}	

	print(paste0("sigma_g :", toString(var_g_est2)))
	#print(paste0("min_rms :", toString(min_rms)))
	#set sigma_g parameter
	parameters[1] <- var_g_est2

	######### estimate sigma_c1 #############
	c1_est2<-maxVar-1/n1-var_g_est2


	if(c1_est2<0){
 		c1_est2<-0.0000000000000000000000000001
 		parameters[1] <- maxVar-1/n1 #c1_est2 sampling error, so set to zero
	}
	#set sigma_c1 parameter
	parameters[2] <- c1_est2
	print(paste0("sigma_c1 ", toString(c1_est2)))

	######### estimate sigma_c2 #############
	max<- - Inf
	c2_est<-0.0000000000000000000000000001
	stats <- cbind(s1, s2)
	#print(dim(stats))
	lik <- c()
	for(i in seq(from=0.0000000000000000000000000001,to=10, by=0.01)){
		temp<-MLE_joint_probability(var_g_est2, c1_est2, i, n1, n2, stats)
		if(temp>max){
			max<-temp
			c2_est<-i
		}
		lik <- c(lik, temp)
	}

	#set sigma_c2 parameter
	parameters[3] <- c2_est
	print(paste0("sigma_c2 ", toString(c2_est)))

	return(parameters)
}


MLE_joint_probability<-function(var_g, var_c1, var_c2, n1, n2, stats){
	cov_matrix=matrix(data=NA, nrow=2, ncol=2)
	cov_matrix[1,1]=var_g+var_c1+1/n1
	cov_matrix[1,2]=var_g
	cov_matrix[2,1]=var_g
	cov_matrix[2,2]=var_g+var_c2+1/n2

	mean_matrix=matrix(data=NA, nrow=2, ncol=1)
	mean_matrix[1,1]=0
	mean_matrix[2,1]=0

	estimate<-0

	prob<-dmvnorm(x=stats, mean = mean_matrix, sigma = cov_matrix, log = TRUE)

	for(i in prob){
 		estimate<-estimate + i
  	}
  	return(estimate)
}







predict_replication <- function(mean, sd, z){
	#calculate predicted replication rate
	lower <-  pnorm(z, mean, sd)
	upper <- 1- pnorm(-z, mean, sd)
	return(lower + upper)
}


#calculate conditional probability of s2 > t | s1 = x under model with no confounding
calc_conditional_no_confounding <- function(s1, n1, n2, sigma_g, z){
	mean <- s1 * (sigma_g)/(sigma_g + 1/n1) 
	sd <- sqrt(1/n2 + sigma_g - (sigma_g*sigma_g)/(sigma_g +1/n1))


#	mean <- (sqrt(n1*n2)*var_g)/(1+n1*var_g)
#	sd <- sqrt(n2*var_g+1-((n1*n2*var_g*var_g)/(n1*var_g+1)))
	lower <-  pnorm(z, mean, sd)
	upper <- 1- pnorm(-z, mean, sd)
	return(lower + upper)
}

#calculate conditional probability of s2 > t | s1 = x under model with confounding
calc_conditional_with_confounding <- function(s1, n1, n2, sigma_g, sigma_c1, sigma_c2, z){
	mean <- s1* ( sigma_g) / (sigma_g + sigma_c1 +1/n1)
	sd <-sqrt( sigma_g + sigma_c2 + 1/n1 - ((sigma_g*sigma_g)/(sigma_g + sigma_c1 +1/n1) ))


#	mean <- (sqrt(n1*n2)*var_g)/(1+n1*var_g+n1*var_c1)
#	sd <- sqrt(n2*var_g+n2*var_c2+1-((n1*n2*var_g^2)/(n1*var_g+n1*var_c1+1)))
	lower <-  pnorm(z, mean, sd)
	upper <- 1-pnorm(-z, mean, sd)
	return(lower + upper)

}

#calculate predicted replication rate under model with no confounding
predict_no_confounding <- function(data, z, sigma_g){
	#currently assume sample size is same for all variants
	n1 <- data[1,3]
	n2 <- data[1,4]
#	sigma_g <- estimate_sigma_g(data[,1], data[,2], n1, n2)
	
	predicted_replication <- 0
	for(i in 1:nrow(data)){
		predicted_replication <- predicted_replication +  calc_conditional_no_confounding(data[i,1], n1, n1, sigma_g, z)

	}

	return(predicted_replication/nrow(data))
}

#calculate predicted replication rate under model with confounding
predict_with_confounding <- function(data, z, sigma_g, sigma_c1, sigma_c2){
	n1 <- data[1,3]
	n2 <- data[1,4]
#	sigma_g <- estimate_sigma_g(data[,1], data[,2], n1, n2)
#	sigma_c1 <- estimate_sigma_c1(data[,1], n1, sigma_g)
#	sigma_c2 <- estimate_sigma_c2(data[,2], n2, sigma_g)

	predicted_replication <- 0
	for(i in 1:nrow(data)){
		 predicted_replication <- predicted_replication +  calc_conditional_with_confounding(data[i,1], n1, n2, sigma_g, sigma_c1, sigma_c2, z)

	}

	return(c(predicted_replication/nrow(data), sigma_g, sigma_c1, sigma_c2))
}

calcReplication <- function(data, z){
	count <- 0
	for(i in 1:nrow(data)){
		if(abs(data[i,2])>abs(z)){
			count  <- count + 1
		}
	}
	return(count/nrow(data))
}



plot_conditional <- function(data, slope, sd, slope_noConfounding, sd_noConfounding, rep){
	maxS <- max(c(abs(data[,1]), abs(data[,2])))*2
	ggplot(data)+
#		scale_color_manual(values=c("grey", "black"))+
		scale_color_manual(values = c("Winner's curse + confounding" = "#ca0020", "Winner's curse"="#0571b0"))+
		geom_point(aes(x=data[,1], y=data[,2], colour=rep))+
		geom_abline(aes(color="Winner's curse"), intercept = 0, slope = slope_noConfounding, color="#0571b0", size=1.25)+
		geom_abline(intercept = -2*sd_noConfounding, slope = slope_noConfounding, linetype=2, color="#0571b0", size=1.25)+
		geom_abline(intercept = 2*sd_noConfounding, slope = slope_noConfounding, linetype=2, color="#0571b0", size=1.25)+
		geom_abline(aes(color="Winner's curse + confounding"), intercept = 0, color="#ca0020", slope = slope, size=1.25)+
		geom_abline(intercept = -2*sd, slope = slope, color="#ca0020", linetype=6, size=1.25)+
		geom_abline(intercept = 2*sd, slope = slope, color="#ca0020", linetype=6, size=1.25)+
		theme_bw()+
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
		ylab("Replication Summary Statistics")+
		xlab("Discovery Summary Statistics")+
		theme(legend.position = 'top',
				text = element_text(size=25),
				plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
		coord_cartesian(xlim = c(-maxS, maxS), ylim = c(-maxS, maxS)) 

	ggsave(paste0(args[2], ".png"))
}

run_mouse_analysis <- function(data){
#	bonferonni <- 0.05/nrow(data)

	if (correction %in% c( "--bonferonni")){
		threshold2 <- threshold/nrow(data)
		z <- qnorm(threshold2/2) 
	}
	else{
		z <- qnorm(threshold/2)
	}

	#variants that are significant in initial replication study
	sig <- subset(data, abs(data[,1])> abs(z))
	#print(dim(sig))
	stopifnot(nrow(sig)>1) #cannot apply method if there is not at least one significant variant

	if (correction %in% c( "--bonferonni")){
		threshold2 <- 0.05/nrow(sig)
		z2 <- qnorm(threshold2/2) 
	}
	else{
		threshold2 <- 0.05/nrow(sig)
		z2 <- qnorm(threshold2/2)
	}
	
	#true replication rate
#	bonferonni2 <- 0.05/nrow(sig)
	r_true <- calcReplication(sig, z2)

	
	#estimate parameters
	s1 <- data[,1]
	s2 <- data[,2]
	n1 <- data[1,3]
	n2 <- data[1,4]
	M <- nrow(data)
	MSig <- nrow(sig)
	Zscore <- abs(z)
	parameters <- estimate_parameters(s1, s2, n1, n2, M, MSig, Zscore)
	sigma_g <- parameters[1]
	sigma_c1 <- parameters[2]
	sigma_c2 <- parameters[3]

	#predicted replication rate with no confounding
	pr_no_confounding <- predict_no_confounding(sig, z2, sigma_g)

	#predicted replication rate with confounding, sigma_g1, sigma_c1, and sigma_c2
	pr_with_confounding <- predict_with_confounding(sig, z2, sigma_g, sigma_c1, sigma_c2)


	#output predicted replication rates and parameters
	write.table(c(r_true, pr_no_confounding, pr_with_confounding), file=args[2], row.names=F, col.names=F, quote=F)

	#make summary plot
#	sigma_g <- pr_with_confounding[2]
#	sigma_c1 <- pr_with_confounding[3]
#	sigma_c2 <- pr_with_confounding[4]
	n1 <- data[1, 3]
	n2 <- data[1, 4]

	var_g <- sigma_g
	var_c1 <- sigma_c1
	var_c2 <- sigma_c2
# 	slope <- (var_g)/(1/n1+var_g+var_c1)
#	sd <- sqrt(var_g+var_c2+1/n2-((var_g^2)/(var_g+var_c1+1/n1)))
#	slope_noConfounding <- (var_g)/(1/n1+var_g)
#	sd_noConfounding <- sqrt(var_g+1/n2-((var_g*var_g)/(var_g+1/n1)))

	slope <- (sqrt(n1*n2)*var_g)/(1+n1*var_g+n1*var_c1)
	sd <- sqrt(n2*var_g+n2*var_c2+1-((n1*n2*var_g^2)/(n1*var_g+n1*var_c1+1)))
	slope_noConfounding <- (sqrt(n1*n2)*var_g)/(1+n1*var_g)
	sd_noConfounding <- sqrt(n2*var_g+1-((n1*n2*var_g*var_g)/(n1*var_g+1)))


    rep <- c()
    for(i in 1:nrow(sig)){
        s2 <- sig[i,2]
        if(abs(s2)>abs(z2)){
            rep <- c(rep, "black")
        }
        else{
            rep <- c(rep , "gray70")
        }
    }
	plot_conditional(sig, slope, sd, slope_noConfounding, sd_noConfounding, rep)

}

run_mouse_analysis(data)
	
