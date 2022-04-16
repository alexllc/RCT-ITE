# code to audit behavior of pseudo R package 

time <- df$os_year
event <- df$death_final
tmax <- NULL
# tmax <- floor(max(time))
	if(is.null(tmax) )
		tmax <- max(time[event==1])

# pseudo mean function
	
	rmtime <- ifelse(time >= tmax ,tmax,time)			#if times are greater than tmax, set them to tmax
	rmdead <- ifelse(time > tmax ,0, event)			#if times are greater than tmax, set censoring status to 0

	if(sum(rmdead)==0)
		stop("no events occured before time 'tmax'")
    
	howmany <- length(rmtime)
    	
    ## preparing the data
  	pseudo <- data.frame(id=1:howmany,time=rmtime,event=rmdead)
  
  	# sort in time, if tied, put events before censoring
  	pseudo <- pseudo[order(pseudo$time,-pseudo$event),]

    	# surv.omit: leave one out
        	howmany <- nrow(pseudo)
	
            td <- pseudo$time[pseudo$event==1]
            lt.temp <- c(td[-1],td[length(td)]+1) # +1 for the last value
            lt <- which(td!=lt.temp)

            Y1 <- matrix(howmany:1,byrow=TRUE,ncol=howmany,nrow=howmany)
            Y2 <- matrix((howmany-1):0,byrow=TRUE,ncol=howmany,nrow=howmany)
            Y <- upper.tri(Y1,diag=FALSE)*Y1+lower.tri(Y2,diag=TRUE)*Y2
            N <- matrix(pseudo$event,byrow=TRUE,ncol=howmany,nrow=howmany)
            Ndiag <- diag(diag(N))
            N <- N - Ndiag
            kmji <- (Y-N)/Y

            km <- t(apply(kmji,1,cumprod))

            tt <- matrix(pseudo$time,byrow=TRUE,nrow=nrow(pseudo),ncol=nrow(pseudo))
            #diag(tt) <- c(diag(tt[-nrow(pseudo),-1]),tmax)
            diag(tt) <- c(0,diag(tt[-1,-nrow(pseudo)]))
            tt <- tt[,pseudo$event==1,drop=FALSE]
            tt <- tt[,lt,drop=FALSE]
            tt <- cbind(rep(0,nrow(pseudo)),tt,rep(tmax,nrow(pseudo)))
            tt <- t(apply(tt,1,diff))

            aje <- which(is.na(km[howmany,]))
            if(length(aje)>0){
                kir <- min(aje)
                km[howmany,kir:ncol(km)] <- km[howmany,kir-1] 
            }

            km <- km[,pseudo$event==1,drop=FALSE]
            km <- km[,lt,drop=FALSE]
            if(!missing(tmax)){
                km <- apply(cbind(rep(1,nrow(pseudo)),km)*tt,1,sum)
            }
            RM.omit <- km

        # surv.tot 
            howmany <- nrow(pseudo)
            
            td <- pseudo$time[pseudo$event==1]
            lt.temp <- c(td[-1],td[length(td)]+1)
            lt <- which(td!=lt.temp)
            
            #km - i
            Y <- howmany:1
            N <- pseudo$event
            
            kmji <- (Y-N)/Y
                
            km <- cumprod(kmji)
            
            if(!missing(tmax)){
                tt <- pseudo$time[pseudo$event==1]
                tt <- tt[lt]
                tt <- c(0,tt,tmax)
                tt <- diff(tt)
            }
            
            #only for deaths, one value per tie
            km <- km[pseudo$event==1]
            km <- km[lt]
            if(!missing(tmax)){
                km <- sum(c(1,km)*tt)
            }
            RM.tot <- km

        pseu <- howmany*RM.tot - (howmany-1)*RM.omit

	
	#RM, all cases
	RM.tot <- surv.tot(pseudo,tmax)

	# pseudo-observations
	pseu <- howmany*RM.tot - (howmany-1)*RM.omit

	#back to original order
	pseu <- pseu[order(pseudo$id)]