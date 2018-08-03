#Compute incidence of cores in a set of events using inclusion metric. If
#the events matrix represents events in a single profile, the result is a row in
#the incidence table.
#
#Input: 
#cores and events are both two-column matrices, with columns named "start" and 
#"end".
#dropevents (character) is a switch for the event-core matching rule. If
#dropevents="Drop" (default), best match is found
#among all pairs of events and cores and the incidence of the best matching core
#is set to that value. Both the event and the core are dropped, and the process
#is repeated until there are no cores left. If dropevents is not "Drop" or 
#"Greedy", only the core is dropped and the process is repeated until there 
#are no cores left. If dropevents="Greedy", the process essentially mimicks the
#greedy implementation of the CORE algorithm, with the order of the cores as 
#specified by the input. But the score for a core is maximized, not summed,
#over all events.
#
#Value: a numerical vector of length equal to nrow(cores), containing matching
#values for the cores.
#
incidence<-function(cores,events,dropevents="Drop",assoc="I"){
	fulltruth<-matrix(ncol=nrow(cores),nrow=nrow(events),data=0)
	if(assoc=="I")for(i in 1:nrow(cores))fulltruth[,i]<-
		(events[,"start"]<=cores[i,"start"]&events[,"end"]>=cores[i,"end"])*
		(cores[i,"end"]-cores[i,"start"]+1)/(events[,"end"]-events[,"start"]+1)
	else if(assoc=="J")for(i in 1:nrow(cores))fulltruth[,i]<-
		pmax(((pmin(events[,"end"],cores[i,"end"])-pmax(events[,"start"],cores[i,"start"])+1)/
		(pmax(events[,"end"],cores[i,"end"])-pmin(events[,"start"],cores[i,"start"])+1)),0)
	else stop("Undefined association")
	coretruth<-rep(0,nrow(cores))
	if(nrow(events)==0)return(coretruth)
	if(dropevents!="Greedy"){
		corelabel<-1:nrow(cores)
		for(i in 1:nrow(cores)){
			bestmatch<-which.max(fulltruth)
			whichcore<-(bestmatch-1)%/%nrow(fulltruth)+1
			whichevent<-(bestmatch-1)%%nrow(fulltruth)+1
			coretruth[corelabel[whichcore]]<-fulltruth[bestmatch]
			if(dropevents=="Drop")fulltruth<-fulltruth[-whichevent,,drop=F]
			fulltruth<-fulltruth[,-whichcore,drop=F]
			corelabel<-corelabel[-whichcore]
		}
	}
	else{ 
		weights<-rep(1,nrow(events))
		for(i in 1:nrow(cores)){
			bestmatch<-which.max(fulltruth[,i]*weights)
			coretruth[i]<-fulltruth[bestmatch,i]*weights[bestmatch]
			weights[bestmatch]<-weights[bestmatch]*(1-fulltruth[bestmatch,i])
		}
	}
	return(coretruth)
}
