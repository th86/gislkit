#Evaluate the extreme value association (EVA) between two vectors
#Reference: http://www.biomedcentral.com/1755-8794/3/51
eva<-function(vec,resp){
	vec.sorted=sort(vec, decreasing=TRUE)
	resp.sorted=resp[names(vec.sorted)]

		eva_vec<-rep(0, length(resp.sorted))
		for(i in 1:length(vec.sorted)){
			x=length(which(resp.sorted[1:i]==0)) #W
			m=length(which(resp.sorted==1)) #W
			n=length(which(resp.sorted==0)) #B
			k=i
			eva_vec[i]=-1*log10( phyper(x, n, m, k, log = FALSE))
			#cat(  phyper(x, n, m, k, log = FALSE), eva[i] , "\n")
		}

		return(eva_vec)
}