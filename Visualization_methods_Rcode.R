
##############################
# Rebecca Worsley Hunt -- 2013.11.20 #
# Wasserman Lab                                #
# University of British Columbia             #
##############################


# source("R_function--CB-plot_TFBS-landscape_TFBS-bi-motif_Dinucleotide-environment.R")

#Examples:
# alias="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/oPOSSUM_results_dinucShuf/Cmyc_hESC/results.txt"
# alias="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/overlap20perc_palind/Nfyb_wgEncodeSydhTfbsHelas3NfybIggrabPk.narrowPeak.1000bp.NFYB-PWM.20percOverlap.top10.seqscores.alt"
# alias="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/overlap20perc_palind/Usf2_wgEncodeSydhTfbsGm12878Usf2IggmusPk.narrowPeak.1000bp.USF1-PWM.20percOverlap.top10.seqscores.alt"
# alias="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/overlap20perc_palind/Srebp1_wgEncodeSydhTfbsHepg2Srebp1InslnStdPk.narrowPeak.1000bp.SREBF1-SREBP1-HEPG2-PWM.20percOverlap.top10.seqscores.alt"
# alias="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Srebp1_example_absScore.alt"
# alias ="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Stat1_wgEncodeSydhTfbsHelas3Stat1Ifng30StdPk.narrowPeak.1000bp.STAT1-PWM.20percOverlap.top10.seqscores.score85.120proxpeakMax.dinuc.align"
# alias ="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Ctcf_sample_1000bp_forAlignment_startOffset.dinuc.align"   # ok
# alias ="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Ctcf_sample_1000bp_forAlignment_mid500.dinuc.align"    # off
# alias ="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Ctcf_sample_1000bp_forAlignment_mid501.dinuc.align"    # ok
# alias ="~/RESEARCH_link/prj_COMPOSITION/RESULTS_2012/scores_and_positions/data_files/Encode_hg19/Ctcf_sample_1000bp_forAlignment_mid501width.dinuc.align"   #no


cat( "Composition-bias plot:", paste("", "CB.plot( data.file, col.X=F, col.Y=F, headerline=T,  col.Y.pval=F, title.plot=\"\", yaxis.lim=F, save.plot.to=F, col.TFnames=F )", sep="\t"),  sep="\n" )
cat(sep="\n")
cat( "TFBS-landscape view:", paste("", "landscape.view( data.file, col.X=F, col.Y=F, score.type=\"relative\", headerline=T,  title.plot=\"\", save.plot.to=F, threshold=F, heuristic.plot=F , resid.fold=2)", sep="\t"), sep="\n" )
cat(sep="\n")
cat( "TFBS-bi-motif view:", paste("", "bimotif.view(  data.file, col.X=F, col.Y=F, headerline=T, title.plot=\"\", save.plot.to=F  )", sep="\t"), sep="\n" )
cat(sep="\n")
cat( "Dinucleotide-envinronment view:", paste("", "dinucleotide.view( data.file, title.plot=\"\", x1.adj=300, x2.adj=300, save.plot.to=F  )", sep="\t"),  sep="\n" )
cat(sep="\n")
cat( "view.help()", sep="\n")
cat(sep="\n")

	
#***************************
#  Composition-bias plot
#***************************
CB.plot = function( data.file, col.X=F, col.Y=F, headerline=T, title.plot="",  yaxis.lim=F, save.plot.to=F, col.Y.pval=F, col.TFnames=F ){
			# require column names except for dinucleotide environment plot
			if( col.X == F | col.Y ==F ){
				print("Column identifiers for col.X and/or col.Y arguments are missing")
				stop();
			}
			print("loading data....")
			#data.file =  read.delim(file=data.file, as.is=T, header=headerline )
	
			# establish x-axis  limits
			if( max(data.file[, col.X], na.rm=T) <=1  ){ 
				xlimit = c( 0, 1)
			}else{
				xlimit = c( 0, 100)
			}
			# determine maximum enrichment score to establish y-axis limits
			if( col.Y.pval == F){
				# substitute "Inf" values that arise 
				max.score = max(data.file[which( is.finite(data.file[, col.Y ]) == T), col.Y ], na.rm=T)  # Inf is a value in oPOSSUM software
				inf.ind = which( is.infinite(data.file[, col.Y ])==T)
				if( length(inf.ind)>0 ){
					if( max.score > 400){ 
							data.file[inf.ind, col.Y ] = max.score +100  
					}else{  
							data.file[inf.ind, col.Y ] = 500  
					}
				}
				if( yaxis.lim==F){
					ylimit = c( min(data.file[, col.Y], na.rm=T), max(data.file[, col.Y], na.rm=T)+0.15*max(data.file[, col.Y], na.rm=T) )
				}else{
					ylimit = c( min(data.file[, col.Y], na.rm=T), yaxis.lim )					
				}
			}else{ # col.Y.pval==T
				if( yaxis.lim==F){
					ylimit = c( 1, 0 )
				}else{
					ylimit = c( 1, yaxis.lim )	
				}
			}
			# PLOTTING
			print("plotting CB plot....")
			if( save.plot.to != F){
				location=paste(save.plot.to, ".pdf", sep="")
				print(paste("saving to... ", location, sep="") )
				pdf( file=location, width=5, height=5, pointsize=12 ) 
			}
			plot(data.file[, col.X ], data.file[, col.Y ], cex=0.5 , cex.axis=1, cex.lab=1, cex.main=0.9, main= title.plot,  xlab="PFM GC content", ylab=paste("enrichment score (column: ", col.Y,")", sep="" ), xlim=xlimit, ylim=ylimit )
			if( col.TFnames != F){
					text( data.file[, col.X ]+xlimit[2]*0.01, data.file[, col.Y ]+ylimit[2]*0.01,  data.file[, col.TFnames] , cex=0.6, srt=0, adj=0)
			}
			if( save.plot.to !=F){
				dev.off()
			}
	}  # end Composition-bias plot
	
#***************************
#  TFBS-landscape view  ... optional heuristic thresholds
#***************************
landscape.view = function( data.file,  col.X=F, col.Y=F, score.type="relative",  headerline=T, title.plot="", save.plot.to=F, threshold=F, heuristic.plot=F , resid.fold=2, score.diff=0.20){
			# require column names except for dinucleotide environment plot
			if( col.X == F | col.Y ==F ){
				print("Columns identifiers for col.X and/or col.Y arguments are missing")
				stop();
			}
			print("loading data....")
			data.file =  read.delim(file=data.file, as.is=T, header=headerline )
			score.type=tolower(score.type)
			if( score.type == "relative"){
				green.score.thr = 85   # the score cutoff for the 2nd histogram that catches high motif scores
				thr.too.high = 95           # the cutoff for motif score thresholds that are accepted
			}else{
				print("PWM scores provided, not relative motif scores")
				mx = max(data.file[, col.Y], na.rm=T)
				green.score.thr = 8  # cant make equivalent across PWMs because don't have full range of motif scores, but hopefully this will do
				thr.too.high =  mx-2  # backing off of the maximum score a bit, hopefully will do
			}
			
			if( threshold == T | heuristic.plot == T){ 
				#--------------------------------------------------------
				# Enrichment zone boundary -- distance from peakMax
				#--------------------------------------------------------
				print("Calculate enrichment boundary.")
				step=5   # bin size
				max.dist = round(max( abs(data.file[, col.X]), na.rm=T ) / step )*step  # round max value up to nearest 5 
				if( max.dist*2 < 600){
					print("less than +/-250bp the peakMax is too short a sequence. There must be more than 50bp of background sequence.")
					stop();
				}
				max.dist = round( (max.dist - max.dist*0.09)/step)*step  # remove a small width from maximum distance  (needed under some circumstances e.g. repeat masking)
				bin.boundary = seq( 0, max.dist-step, by=step)  # bin value will become the upper limit of the bin
				count.bin = numeric(0) ;
				countfreq = numeric(0) ;
				for(  x in bin.boundary ){
					count.bin = c(count.bin,  length( subset(data.file , abs(data.file[,col.X]) >= x & abs(data.file[,col.X]) < (x+step) )[, 1]) )  
				}
				countfreq = count.bin/length(data.file[,1])    # fraction of seq in a bin relative to all seq
				bin.boundary = c(bin.boundary[-1], max.dist )	# remove index 1, value 0, and add max.dist as the last value ( e.g. max.dist=475,  bin is 470-475bp)
				names( bin.boundary) = seq( 0, max.dist-step, by=step)  # index name is the lower limit of the boundary

				res.fold = resid.fold  # the "residual distance fold adjustment" for the regression line (shift along the y-axis)
				distal.range=c(200, max.dist)  
				left.distal.bin.index = floor(distal.range[1]/step)  
				right.distal.bin.index = length(bin.boundary)
				x.bins.distal= bin.boundary[ left.distal.bin.index : right.distal.bin.index]
				y.fractionSeqs.distal = countfreq[ left.distal.bin.index : right.distal.bin.index]
				count.line.model = lm( y.fractionSeqs.distal ~ x.bins.distal)  # model on the distal bins
				count.coef.line = coef( count.line.model)  # the  y intercept and slope
				max.abs.res = max( count.line.model$residuals[order(count.line.model$residuals)][2:(length(count.line.model$residuals)-2)]   ) 	# 3rd largest residual
				# generate y-values of the line from the regression model (x-value is the bins)
				dist.bin = distal.range[1]/step		
				line.y.values = count.coef.line[2] * bin.boundary + count.coef.line[1]
				adjusted.y.values = line.y.values + (max.abs.res*res.fold )  # adjust regression line y-values be 2fold the 3rd residual
				
				# Decision on where to place enrichment threshold. Look at how consistently values are above the adjusted regression line i.e. "allowance line"
				flag=1
				bin.index.above2fold = which( countfreq[1: (dist.bin )] >=  adjusted.y.values[1:dist.bin]  )  # bins above or equal to the adjusted regression line
				if(length(bin.index.above2fold ) == 0){
					bin.index.above2fold==0   # this will result in one bin boundary.... although maybe I should indicate when zero
					print("WARNING: poor signal, boundary can't be set")
					flag=0
				}else{ 
					bin.index.below2fold = which( countfreq[1:(dist.bin )] <=  adjusted.y.values[1:dist.bin]   ) # bins BELOW OR EQUAL 2fold adjusted regression
					# Sometimes the set of bins closest to the background regions, goes varies above-below-above the adjusted regression line. 
					# The below decides if there are enough bins above the adjusted line to set the enrichment threshold, or if we have to drop some of the bins before setting the threshold 
					for( ind in length(bin.index.above2fold):1 ){ 
						flag=1
						last.bin.above = bin.index.above2fold[ind]
						uncertain.bins.below = bin.index.below2fold[ which(bin.index.below2fold < last.bin.above) ]   # the bins below the line that are intermixed with bins above line
						count.mixed.bins = last.bin.above-uncertain.bins.below[1] +1 # number of bins that are mixture of up and down
						# if too great a proportion of mixed bins are below adjusted2fold line, then the bin.index.above2fold needs to be adjusted to only the bins where everything beyond is above line. 
						if( length(uncertain.bins.below) > count.mixed.bins*0.4  && !is.na(count.mixed.bins)  &&  tail(bin.index.above2fold, n=1) > 2 ){ 
							 flag = 0;
							 if( ind==1 && uncertain.bins.below[1] == 1){  # if you are on last point, then all is lost, nothing enriched
							 	bin.index.above2fold=0 
							 	print("WARNING indecisive signal: boundary can't be set")
							 	break;
							}else{
						 		 bin.index.above2fold = bin.index.above2fold[1:(ind-1)]
						 	}
						 	ind=length(bin.index.above2fold)
						}
						if( flag==1){   # if came through for loop with flying colours, all is well, stop here
							break; 
						}
				  	}
				}
				enrichment.threshold = step*tail(bin.index.above2fold, n=1) 
				if( enrichment.threshold < 7){   # a motif width of <14
					enrichment.threshold = NA 
				}
				if(flag == 0){ 
					enrichment.threshold = NA 
				}
				#----------------------------------------------------------
				# Calculate motif score threshold. 
				#----------------------------------------------------------
			 	print("Calculate lower motifscore threshold")
				if( ! is.na(enrichment.threshold) ){
					#  collect a background range equal to the enrichment zone width e.g. 150bp of foreground, thus collect 150bp of background to analyze
					# pull the background from multiple samplings (i.e. windows) of the distal zone rather than using one window because distal zone can be variable
					# then evaluate the motifscores in the enrichment zone vs the background control region
					bckgrnd.buffer = 50  # background starts 50bp outside of enrichment zone
					win.width = 5             # 5bp per window, will combine multiple windows to make background
					distalrange= max.dist - (enrichment.threshold+ bckgrnd.buffer)  
					numwin= enrichment.threshold/win.width
					lastbit = win.width*(numwin%%1)  # instances where the number of windows is not a whole number due to enrichment.threshold not evenly divisable by 5
					numwin=floor(numwin)  # number of background windows
					starts=seq(enrichment.threshold+bckgrnd.buffer, max.dist, by=floor(distalrange/numwin) )  # the distance that each 5bp window will start from
					distal=matrix(nrow=0, ncol=0)
					for( ind in 1:numwin ){
						distal=rbind(distal, data.file[which(abs(data.file[, col.X]) >= starts[ind] & abs(data.file[, col.X]) < starts[ind]+5 ), ]  )
					}
					if( lastbit > 0){  # get the last few basepairs that aren't quite a full 5bp
						distal=rbind(distal, data.file[which(abs(data.file[, col.X]) >= (starts[length(starts)]+5*2) & abs(data.file[, col.X]) < (starts[length(starts)]+5*2+lastbit) ), ]  )
					}
					if( score.type=="relative"){
						bin.step=1   # one relative motif score point is the bin size e.g. bin 85, bin 86, bin 8	
					}
					if( score.type == "pwm"){
						bin.step=0.5   # one relative motif score point is the bin size e.g. bin 85, bin 86, bin 87						
					}
					min.score=floor(min( data.file[, col.Y ]))-bin.step 
					max.score=ceiling(max( data.file[, col.Y ])) 
					range.scores = seq(min.score, max.score, by= bin.step) #the upper limits of the bins	
					# get bin counts in steps of 1 motif score for central seqs (foreground) and distal seqs (background)
					central = data.file[which( abs(data.file[, col.X ]) <= enrichment.threshold ), ]
					countcentral=numeric(0)
					countdistal=numeric(0)
					for( score in range.scores ){ 
						countcentral = c(countcentral, length(which(central[, col.Y ] > score & central[, col.Y ] <= score+1 )) ) 
						countdistal = c(countdistal, length(which(distal[, col.Y] > score & distal[, col.Y ] <= score+1 )) )
					}	
									
					names(countcentral)= range.scores			
					names(countdistal)= range.scores
					# put the max distal value from the set of control bins to right of (and including) the bin of interest, into elements of array corresponding to bin of interest.
					# If walk along the plot of bin counts, can see that once past the global maximum, the max distal value IS the bin of interest
					# Result is that low score bins will be more stringently assessed than high score bins. 
					maxfwdcountdistal = rep.int( 0, length(range.scores) )
					for( ind in 1:length(countdistal) ){
						maxfwdcountdistal[ind] = max(countdistal[ind:length(countdistal)], na.rm=T )
					}
					names(maxfwdcountdistal) = range.scores
					frac.test = score.diff   #  require the central count to exceed the control count by 20% of the global maximum of the set of controls with bin names greater than the bin of interest 
					fraction = (countcentral - countdistal)/maxfwdcountdistal
					flag="beginthr"
					frac.flag=numeric(0)
					for(index in 1:length(fraction) ){
						if( is.na(fraction[index]) | fraction[index]<= frac.test ){ 
							frac.flag=c(frac.flag, -1)   # bin of interest is below distal count either 2 bins back or 2 bins forward
						}else{
							frac.flag=c(frac.flag,0)
						}
						names(frac.flag)[length(frac.flag)] = names(fraction)[index]
					}	
					#summarize 5consecutive bins around bin of interest. 0 = ok, <0 (e.g. -1,-2,-3,-4) will indicate how many failures of central to exceed distal
					sumry=numeric(0)
					for(ind in 3:(length(frac.flag)-2) ){
						sumry=c(sumry, frac.flag[ind-2]+frac.flag[ind-1]+frac.flag[ind]+frac.flag[ind+1]+frac.flag[ind+2] )
						names(sumry)[length(sumry)] = names(frac.flag)[ind]
					}
					ind.thr = names(sumry)[tail( which(sumry ==0), n=1) ] 
					# if there is not a value of 0 found in sumry, then look for a value of -1 (only one failure)
					if(length(ind.thr)==0){
						ind.thr = names(sumry)[tail( which(sumry ==-1), n=1) ] 
					}
					# if there is not a value of -1 found in sumry, then flag "fail" (there were 2 or more failures)
					if(length(ind.thr)==0){ # if no 0 or -1, then set ind.thr at right
						flag="failthr"
						thr.motifscore=NA
					}
					if( flag != "failthr"){
						while(sumry[ind.thr] > -2){ 
							ind.thr=as.character(as.numeric(ind.thr) -bin.step); #walk down until reach -2
							for(ind in bin.step*2:bin.step){
								if( sumry[as.character(as.numeric(ind.thr)-ind)]== -1 ){ #if the next element below -2 is a -1, then keep it and keep walking (rare case). -2 is ok on its own, avoid a pair or more of them
									ind.thr=as.character(as.numeric(ind.thr)-ind); 
								}
							}
						}
						ind.thr=as.character(as.numeric(ind.thr)) 
						# make sure the chosen threshold isn't a choice where central is less than control, if it is, then go one step forward to next bin 
						while(frac.flag[ind.thr] ==-1){ 
							ind.thr=as.character(as.numeric(ind.thr)+bin.step);
						}
						thr.motifscore = as.numeric(ind.thr)
						if( thr.motifscore > thr.too.high ){  #if it gets this high (relative score 95), then its because there is nothing to threshold
							thr.motifscore=NA
						}
					} # end of test: flag not failthr			
				}else{  #boundary is 0
					thr.motifscore = NA
				}
				results=  c( enrichment.threshold , enrichment.threshold*2, thr.motifscore ) 
				names(results) = c("absolute_enrichment_boundary", "enrichment_zone_width", "motifscore_threshold")
			}	# end thresholds
			#----------------------------------------------------------
			#  Plotting landscape view (and thresholds lines)
			#----------------------------------------------------------
			print("plotting TFBS-landscape view....")
			if( save.plot.to !=F){
				location=paste(save.plot.to, ".pdf", sep="")
				print(paste("saving to... ", location, sep="") )
				pdf(file=location, width=7.5, height=4.5, pointsize=12 ) 
			}
			par(mfrow=c(1,2) )	
			# set x-axis limits
			xlimit = c( min(data.file[, col.X ], na.rm=T), max(data.file[, col.X ], na.rm=T) ) 
			# set y-axis limits for plot 1
			ylimit.score=c(min( data.file[, col.Y ]), max( data.file[, col.Y ]))  
			par(mar=c(4.5, 4.5, 2.5, 0.5))  
			plot( data.file[, col.X ], data.file[, col.Y ], cex=0.1, cex.main=0.7, main=title.plot, xlab="motif distance to peakMax", ylab="motif score", xlim=xlimit, ylim=ylimit.score ) 
			if( threshold == T){
				abline(v=c(-1*enrichment.threshold, enrichment.threshold), lty=2, col="blue", lwd=1.5)
				abline(h=thr.motifscore, lwd=1.5, lty=2, col="blue")
			}
			# green line, histogram of sequences with motif relative score above 85
			par(mar=c(4.5, 4.5, 2.5, 2.5))  
			brknum=(round((max(data.file[,"distance1"])-min(data.file[,"distance1"]))/50)*50) / 5    # 5bp resolution 
			high.score.data = data.file[which(data.file[, col.Y]>= green.score.thr ), ]
			high.score =hist( plot=F, high.score.data[, col.X ], breaks= brknum )
			max.dens = max(high.score$density, na.rm=T) 
			if( max.dens <= 0.006){
				ylimit.green = c(0, 0.006 )
			}else{
				ylimit.green = c(0, max.dens )
			}
			plot(high.score$mids, high.score$density, type="l",   cex=0.1, cex.main=0.7, main= title.plot, xlab="", ylab="", xlim=xlimit, ylim= ylimit.green , col="green",  yaxt="n"); 
			axis(side=4, col="darkgreen"  )
			# black line, histogram of all sequences
			par(new=T)  
			brknum=(round((max(data.file[,"distance1"])-min(data.file[,"distance1"]))/50)*50) / 2   # 2bp resolution
			all.score = hist( plot=F, data.file[, col.X ], breaks= brknum ) 
			max.dens = max( all.score$density, na.rm=T) 
			if( max.dens <= 0.006){
				ylimit = c(0, 0.006 )
			}else{
				ylimit = c(0, max.dens )
			}	
			plot(all.score$mids, all.score$density, type="l", cex=0.1, cex.main=0.7, main="", xlab="motif distance to peakMax", ylab="probability density",  xlim=xlimit,  ylim=ylimit, col="black" )
			legend("topleft", legend=c("all motifs", paste("high scoring (>", round(green.score.thr, digits=1), ")", " motifs" ,sep="") ), col=c("black","green"), cex=0.7, lty=1,lwd=1.5, bty="n")
			if( threshold == T){
				abline(v=c(-1*enrichment.threshold, enrichment.threshold), lty=3, col="blue", lwd=1.5)
			}
			if( save.plot.to !=F){
				dev.off()	
			}
			#----------------------------------------------------------
			#  Plotting heuristic decision plots
			#----------------------------------------------------------
			if( heuristic.plot == T ){
				print("plotting heuristic decisions....")
				if( save.plot.to != F){
					location=paste(save.plot.to, ".heuristicDecision.pdf", sep="")
					print(paste("saving to... ", location, sep="") )
					pdf(file=location, width=11.25, height=4.5, pointsize=12 )  	
				}else{
					devAskNewPage(ask=T) 
				} 
				par(mfrow=c(1,3))
				par(mar=c(5.7,5,3.5,2) )
				par( oma = c( 0, 0, 1, 0 ) ) #side 3 has space for title
				# set x-axis limits;  distance from peakMax
				xlimit = c( min(data.file[, col.X ], na.rm=T), max(data.file[, col.X ], na.rm=T) ) 
				# set y-axis limits for plot 1
				ylimit.score=c(min( data.file[, col.Y ]), max( data.file[, col.Y ])) 
				# plot 1
				plot( data.file[, col.X ], data.file[, col.Y ], cex=0.1, cex.main=1.2, main=title.plot, cex.lab=1.2, cex.axis=1.2, xlab="motif distance to peakMax", ylab="motif scores", xlim=xlimit, ylim=ylimit.score )
				abline(v=c(-enrichment.threshold, enrichment.threshold), lwd=1.5, lty=1, col="blue")
				abline(h=thr.motifscore, lwd=1.5, lty=1, col="blue")
				# plot 2
				ymax=0.06
				if(max(countfreq, na.rm=T)>0.06){
					ymax=max(countfreq[1:10],na.rm=T)
				}
				plot( bin.boundary, countfreq , pch=16, cex=0.5, main="",  cex.lab=1.2, cex.axis=1.2, xlab="", ylab="proportion of peaks in bin",  ylim=c(0, ymax), xaxt="n")
				title(main="Decision: enrichment zone absolute boundary", cex.main=1 )
				title(xlab="upper limit of binned absolute \ndistance to peakMax (bin = 5bp)" , cex.lab=1.2, cex.axis=1.2, line=4.5)
				b.labels=bin.boundary[seq(1,length(bin.boundary), by=5) ];  
				axis(side=1, labels= b.labels, at=b.labels, las=2, cex.axis=1.2, cex.lab=1.2 )
				abline( v =distal.range[1] ,  lty=1, col="black",  lwd=1.3)
				abline(count.coef.line , lty=1, col="green",  lwd=1.1)
				abline(c(count.coef.line[1]+max.abs.res*res.fold, count.coef.line[2]), lty=2, col="dodgerblue" ,  lwd=1.1) 
				abline( v =enrichment.threshold ,  lty=1, col="blue",  lwd=1.3)
				legend("topright", legend=c( "enrichment boundary", "allowance line", "distal values to right of line", "regression of distal values" ), col=c("blue", "dodgerblue", "black", "green"), lty=c(1,2,1,1), cex=1, bty="n" )
				#plot 3
				if( !is.na(enrichment.threshold) ){
					xlimit=c( as.numeric(head(names(countcentral), n=1)) , as.numeric(tail(names(countcentral), n=1)) )
					ylimit=c(0, max(countcentral) )
					plot(names(countcentral), countcentral, cex=0.5, pch=16, main="", cex.lab=1.2, cex.axis=1.2, xlab="upper limit motif score of bin", ylab="frequency of peaks",  type="b" , xlim=xlimit, ylim=ylimit )
					par(new=T); 
					plot(names(countdistal), countdistal, cex=0.5, pch=16, main="",  cex.lab=1.2, cex.axis=1.2,  type="b", xlab="", ylab="", xlim=xlimit, ylim=ylimit, col="red",  xaxt="n", yaxt="n")
					abline(v= thr.motifscore , lty=1 , col="blue" )
					legend("topleft", legend=c("enrichment zone", "control"), text.col=c("black","red"), col=c("black","red"), cex=1, pch=16, bty="n")
					title(main="Decision: motif score threshold", cex.main=1 )
				}
				if( save.plot.to != F){
					dev.off()
				}else{
					devAskNewPage(ask=F) 
				}
			}  # end of decision to plot heuristic decision plots
			if( threshold==T | heuristic.plot==T){
				return(results)
			}
	}  # end TFBS-landscape view 
	
#***************************
#  TFBS-bi-motif view
#***************************
bimotif.view =function( data.file,  col.X=F, col.Y=F, headerline=T, title.plot="", save.plot.to=F ){
			# require column names except for dinucleotide environment plot
			if( col.X == F | col.Y ==F ){
				print("Columns identifiers for col.X and/or col.Y arguments are missing")
				stop();
			}
			print("loading data....")
			data.file =  read.delim(file=data.file, as.is=T, header=headerline )
			print("plotting TFBS-bi-motif view....")
			if( save.plot.to !=F){
				location=paste(save.plot.to, ".pdf", sep="")
				print(paste("saving to... ", location, sep="") )
				pdf(file=location, width=8.5, height=4, pointsize=12 ) 
			} 
			par(mfrow=c(1,2) )
			resolution=5
			min.x =round(min(data.file[, col.X ] - data.file[, col.Y ], na.rm=T)/100)*100
			max.x =round(max(data.file[, col.X ] - data.file[, col.Y ], na.rm=T)/100)*100
			min.y =round(min(data.file[, col.Y ], na.rm=T)/100)*100
			max.y =round(max(data.file[, col.Y ], na.rm=T)/100)*100
			xlimit=c( min.x, max.x)
			ylimit=c( min.y , max.y )
			plot( data.file[, col.Y ] - data.file[, col.X ], data.file[, col.X ], cex=0.1, cex.axis= 1, cex.lab= 1, cex.main=0.7, main= title.plot, xlab="distance of motif2 from motif1", ylab="motif1 distance to peakMax" , xlim=xlimit, ylim=ylimit ) 
			hist.values= hist( plot=F, data.file[, col.Y ] - data.file[, col.X ],  breaks= max.x/resolution )
			if( max(hist.values$density, na.rm=T) > 0.006){
				ylimit.hist=c(0, max(hist.values$density, na.rm=T) )
			}else{
				ylimit.hist=c(0, 0.006)
			}
			hist( data.file[, col.X ] - data.file[, col.Y ],  breaks=1000/resolution, freq=F, cex.axis= 1,  cex.lab= 1, col="gray90", cex.main=0.8, main="", xlab="distances between motifs", xlim=xlimit, ylim=ylimit.hist )
			if( save.plot.to !=F){
				dev.off()	
			}	
	}	# end TFBS-bi-motif view
	
#***************************
#  Dinucleotide-environment view
#*************************** 
dinucleotide.view = function( data.file,  headerline=T, title.plot="", x1.adj=300, x2.adj=300, save.plot.to=F ){
			print("loading data....")
			data.file =  read.delim(file=data.file, as.is=T, header=T, sep="," , row.names=1 )
			if( save.plot.to != F){
				location=paste(save.plot.to, ".pdf", sep="")
				print(paste("saving to... ", location, sep="") )
				pdf(file=location, width=5.5, height=7, pointsize=12 ) 
			}
			motif.loc=numeric(0)
			# get approx. motif location (using maximum of dinucleotide) to organize x-axis labels
			for(ind in 1:16){
				len = length(data.file[1,])
				span = seq( len*0.10, len-len*0.10, by=1)   # use span to remove edges
				max.val = max(data.file[ind+1, span ]/data.file[1, span])
				max.loc = which( data.file[ind+1, ]/data.file[1,] == max.val )
				if(length(max.val)>0){
					if(max.val > 0.4 ){
						motif.loc=c(motif.loc, max.loc)
					}
				}else{
					motif.loc=len/2
					print("no dinucleotide was enriched in >40% of the sequences")
				}
			}
			mid=round(mean(motif.loc, na.rm=T ))
			colours=c("aquamarine", "bisque4", "blue", "blueviolet", "red", "chartreuse3", "chocolate", "plum", "darkgoldenrod2", "deeppink3", "dodgerblue1", "forestgreen", "indianred1", "lightblue2", "cyan4","black")  # tan3  yellowgreen  greenyellow
			names(colours )=c("AA","GG", "CC", "TT", "AG",  "AC", "AT*", "GA", "GC*", "GT", "CA", "CG*", "CT", "TA*", "TG", "TC" )
			xlimit=c(mid-x1.adj, mid+x2.adj)
			ylimit=c(0, 1);
			for( ind in c(1:16) ) {   # 16 dinculeotides = 16 plots
				plot( 1:length(data.file[1,] ), data.file[ ind +1, ]/data.file[1, ],  cex.main=1.0, main="", xlab=xlabel, ylab=ylabel, cex=0.4,  col=colours[ ind ], xlim=xlimit, ylim= ylimit, type="b" , xaxt="n", yaxt="n")   
				if(ind < 16){ par(new=T ) } 
			}
			title(main=title.plot, xlab="nucleotide position", ylab="dinucleotide proportion",  line=2.5 )
			axis(side=1, at=seq(xlimit[1], xlimit[2], by=round((x1.adj+x2.adj)*0.10) ), labels=seq(xlimit[1]-mid, xlimit[2]-mid, by=round((x1.adj+x2.adj)*0.10)) )
			axis(side=2)
			legend.colours =colours[c("AA","TT", "GG","CC","AG", "CT", "AC", "GT",   "GA","TC",    "CA","TG",   "AT*","TA*", "CG*", "GC*" )]
			legend("topleft", legend=names(legend.colours), col= legend.colours, bty="n", pch=15, cex=0.9)
			if( save.plot.to !=F){
				dev.off()	
			}	
	}  # end Dinucleotide-environment view	



#-----------------------------------
# Help Function
#-----------------------------------
view.help = function(){
	options(width=120)
	cat("Usage: ", sep="\n")
	cat( "Note:  non-logical and non-numeric arguments should be in double quotes e.g. title.plot=\"My plot 1\"", sep="\n" )
	cat( "       logical arguments are T for true and F for false", sep="\n" )
	cat("", sep="\n")
#***************************
#  Composition-bias plot
#***************************
	cat("", sep="\n")
	cat( "-------------------------", "Composition-bias plot", "-------------------------", "CB.plot(  data.file, col.X=F, col.Y=F, headerline=T,  col.Y.pval=F, title.plot=\"\", yaxis.lim=F, save.plot.to=F, col.TFnames=F )",  sep="\n" )
	cat("", sep="\n")
	cat( "Arguments:", sep="\n")
	cat( "data.file"); writeLines(strwrap("file location (absolute path and file name)",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "col.X");  writeLines(strwrap("column number or name of PFM GC composition",  indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t") )
	cat( "col.Y");  writeLines(strwrap("column number or name of enrichment scores",  indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t") )
	cat( "headerline"); writeLines(strwrap("logical (T or F), does your file have column names?  (default=T)",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "col.Y.pval");  writeLines(strwrap("logical (T or F), are enrichment scores p-values?  (default=F)",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "title.plot");  writeLines(strwrap("set the main title of the plot (default=\"\")",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "yaxis.lim");  writeLines(strwrap("the upper limit of the y-axis (default=F)",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "save.plot.to");  writeLines(strwrap("set the location (absolute path) and filename  e.g. /dir1/dir2/my_plot_name  to save the plot as a PDF file at the given location  (default=F, does not save plot)",  indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat( "col.TFnames"); writeLines(strwrap("the column number or name for TF names, if they are to be printed on the plot (default=F)",  indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat("", sep="\n")

#***************************
#  TFBS-landscape view  ... optional heuristic thresholds
#***************************
	cat("-------------------------", "TFBS-landscape", "-------------------------",  "landscape.view( data.file, col.X=F, col.Y=F, score.type=\"relative\", headerline=T, title.plot=\"\", save.plot.to=F, threshold=F, heuristic.plot=F, resid.fold=2)", sep="\n" )
	cat("", sep="\n")
	cat( "Arguments:", sep="\n")
	cat( "data.file");  writeLines(strwrap("file location (absolute path and file name)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "col.X");   writeLines(strwrap("column number or name of distances between motif to peakMax", indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t") )
	cat( "col.Y");  writeLines(strwrap("column number or name of motif PWM scores (not p-values)", indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t") )
	cat( "score.type");   writeLines(strwrap("type of motif score: \"relative\" (default) or \"pwm\" (default=\"relative\"). Motif relative scores range from 0-100 (they are derived from motif PWM scores); motif PWM scores vary with the PWM, but might range -50 to +20 ", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "headerline");  writeLines(strwrap("logical (T or F), does your file have column names?  (default = T)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "title.plot");   writeLines(strwrap("set the main title of the plot (default=\"\")", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t") )
	cat( "save.plot.to");   writeLines(strwrap("set the location (absolute path) and filename  e.g. /dir1/dir2/my_plot_name  to save the plot as a PDF file at the given location  (default=F, does not save plot)", indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat( "threshold");   writeLines(strwrap("logical (T or F), option to calculate and indicate the enrichment boundary on the plot (default=F). Return value is a vector: c(enrichment zone boundary, enrichment zone width, PWM motif score threshold)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "heuristic.plot");  writeLines(strwrap("logical (T or F), option to calculate and plot the bins used in the heuristic decision (default=F). Return value is a vector: c(enrichment zone boundary, enrichment zone width, PWM motif score threshold)", indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat( "resid.fold"); writeLines(strwrap( "(default = 2 ), pertains to the enrichment zone boundary decision: the adjustment of the allowance line from the regression line is 2-fold (default) the 3rd largest residual value, the allowance line is used to separate the zone of motif enrichment from the background trend (set parameter: heuristic.plot=T, to view the position of the allowance line)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat("score.diff");  writeLines(strwrap("(default = 0.20,  i.e. 20%), pertains to the motif score threshold decision: the required minimal difference between the count of scores in the enrichment zone versus the count of scores in the control region (set parameter: heuristic.plot=T, to view the enrichment zone counts versus the control region counts)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat("", sep="\n")

#***************************
#  TFBS-bi-motif view
#***************************
	cat("-------------------------",  "TFBS-bi-motif", "-------------------------", "bimotif.view( data.file, col.X=F, col.Y=F, headerline=T, title.plot=\"\", save.plot.to=F)", sep="\n" )
	cat("", sep="\n")
	cat( "Arguments:", sep="\n")
	cat( "data.file"); writeLines(strwrap("file location (absolute path and file name)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "col.X");  writeLines(strwrap("column number or name for distances between motif1 and the peakMax", indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t\t") )
	cat( "col.Y");  writeLines(strwrap("column number or name for distances between motif2 and the peakMax", indent=0, exdent=0, initial="\t\t\t", prefix="\t\t\t\t") )
	cat( "headerline"); writeLines(strwrap("logical (T or F), does your file have column names?  (default=T)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "title.plot");  writeLines(strwrap("set the main title of the plot (default=\"\")", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "save.plot.to");  writeLines(strwrap("set the location (absolute path) and filename  e.g. /dir1/dir2/my_plot_name  to save the plot as a PDF file at the given location  (default=F, does not save plot)", indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat("", sep="\n")

#***************************
#  Dinucleotide-environment view
#***************************	
	cat("-------------------------",  "Dinucleotide-envinronment", "-------------------------", "dinucleotide.view( data.file, title.plot=\"\", x1.adj=300, x2.adj=300, save.plot.to=F)",  sep="\n" )
	cat("", sep="\n")
	cat( "Arguments:", sep="\n")
	cat( "data.file"); writeLines(strwrap("file location (absolute path and file name)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "headerline"); writeLines(strwrap("logical (T or F), does your file have column names?  (default=T)", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "title.plot");  writeLines(strwrap("set the main title of the plot (default=\"\")", indent=0, exdent=0, initial="\t\t", prefix="\t\t\t\t") )
	cat( "x1.adj, x2.adj");  writeLines(strwrap("the +/- distance from the aligned motif to display in the plot (default x1.adj=x2.adj=300,  i.e. 600bp x-axis)", indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat( "save.plot.to");  writeLines(strwrap("set the location (absolute path) and filename  e.g. /dir1/dir2/my_plot_name  to save the plot as a PDF file at the given location  (default=F, does not save plot)", indent=0, exdent=0, initial="\t", prefix="\t\t\t\t") )
	cat("", sep="\n")	
}

