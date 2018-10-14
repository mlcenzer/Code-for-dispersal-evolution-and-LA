setwd('~/Dropbox/sbb_dispersal/displa')
rm(list=ls())
source('src/initialize.R')

#make test landscape
#prms.ls <- make.prms.ls(ltc.size=128,
#                        acl=c(0.001,0.1),
#                        frac.suitable=c(0.3,0.3))
 
#load prms.ls.list
make_edge_list<-function(save.path, dir){
	load(file.path(save.path, 'prms.ls.list.RData'))

	summarize_edges<-function(prms.ls){                    
	#Separate into landscape types
	H1o<-prms.ls$landscape[,,1]
	H2o<-prms.ls$landscape[,,2]
	H3<-H1o*H2o
	H1<-H1o-H3
	H2<-H2o-H3
	H0<-abs((H1+H2+H3)-1)

	#note that every upper shared edge is also a lower shared edge for the other cell. Within a cell type, careful not to double-count edges.

	calc_edges_within<-function(H){
		up_down_edges<-sum(H[c(prms.ls$ltc.size,1:(prms.ls$ltc.size-1)),]*H)
		right_left_edges<-sum(H[,c(prms.ls$ltc.size,1:(prms.ls$ltc.size-1))]*H)
		edges<-sum(up_down_edges+right_left_edges)
		edges
	}

	#If run within habitat type returns double the edges
	calc_edges_between<-function(HA,HB){
		up_edges<-sum(HA[c(prms.ls$ltc.size,1:(prms.ls$ltc.size-1)),]*HB)
		down_edges<-sum(HA[c(2:prms.ls$ltc.size, 1),]*HB)
		right_edges<-sum(HA[,c(prms.ls$ltc.size,1:(prms.ls$ltc.size-1))]*HB)
		left_edges<-sum(HA[,c(2:prms.ls$ltc.size,1)]*HB)
		edges<-sum(up_edges + down_edges + right_edges + left_edges)
	edges
	}

	#null edges of each type should be proportional to the abundance of each habitat type
	#Note that within a habitat type this is also going to double-count edges
	null_edges<-function(HA, HB){
		edge_total<-sum(HA)*4
		prop_HB<-sum(HB)/(((prms.ls$ltc.size)^2)-1)
		exp_edges<-edge_total*prop_HB
	exp_edges
	}

#add calculation for with H1 to without H1, with H2 to without H2
#from shared patch to a patch that doesn't include your resource?

	all_habitat<-abind(H0,H1,H2,H3, (H1+H3), (H2+H3), (H0+H2), (H0+H1), (H1+H2+H3), along=3)

	edge_summary<-apply(all_habitat, 3, function(y) apply(all_habitat, 3, function(x) (calc_edges_between(x, y)-null_edges(x,y))/(sum(x)*4))) #consider subtracting instead of dividing?

	colnames(edge_summary)<-c('H0', 'H1', 'H2', "H3", "with_H1", "with_H2", "without_H1", "without_H2", "all_H")
	rownames(edge_summary)<-c('H0', 'H1', 'H2', "H3", "with_H1", "with_H2", "without_H1", "without_H2", "all_H")
edge_summary
}

	edge_list<-lapply(prms.ls.list, summarize_edges)

	save(edge_list, file=sprintf('~/Dropbox/sbb_dispersal/displa/%s/%s.RData', save.path, dir))
}

make_edge_list(save.path<-'saved/long runs/h1.acl-0.01;h1.frac.s-0.3', dir='edge_list')


#normalize by maximum possible number of edges given the number of cells of each type? Probably not, but I had fun solving it so I'm leaving it here for now.

#max_edges_within<-function(H){
#	cells<-sum(H)
#	base<-floor(sqrt(cells))
#	base_edges<-2*base*(base-1)
#	remainder<-cells-base*base
#	if(remainder>0){
#		if(remainder<=base) add_edges<-remainder*2-1
#		if(remainder>base) add_edges<-remainder*2-2
#		edges<-base_edges+add_edges
#	}
#	else edges<-base_edges
#edges
#}	

#max_edges_within(H3)

