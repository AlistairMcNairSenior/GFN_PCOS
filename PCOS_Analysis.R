
# Clear any old objects
rm(list=ls())

# Load the libraries
library(mgcv)
library(sp)
library(fields)
library(gplots)
library(mixexp)

# A function created to find the outer perimeter over which the surface should be fitted
findConvex<-function(x, y, rgnames, res=101, x.limits=NA, y.limits=NA){
	hull<-cbind(x,y)[chull(cbind(x,y)),]	
	# Either use the specifiec limits for the grid, or use pretty
	if(length(which(is.na(x.limits) == T)) > 1){
		px<-pretty(x)	
	}else{
		px<-x.limits	
	}
	if((length(which(is.na(y.limits) == T)) > 1)){
		py<-pretty(y)	
	}else{
		py<-y.limits	
	}	
	# Get the matrix
	x.new<-seq(min(px, na.rm=T),max(px, na.rm=T),len=res)
	y.new<-seq(min(py, na.rm=T),max(py, na.rm=T),len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}



# A function created to find the outer perimeter over which the surface should be fitted for proportional data
findConvex.prop<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	x.new<-seq(0,1,len=res)
	y.new<-seq(0,1,len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}



fields.colors<-function (n, alpha = 1) 
{
    if ((n <- as.integer(n[1L])) > 0) {
        j <- n%/%3
        k <- n%/%3
        i <- n - j - k
        c(if (i > 0) hsv(h = seq.int(from = 43/60, to = 30/60, 
            length.out = i), alpha = alpha), if (j > 0) hsv(h = seq.int(from = 12/60, 
            to = 9/60, length.out = j), alpha = alpha), if (k > 
            0) hsv(h = seq.int(from = 5/60, to = 0, length.out = k), 
            alpha = alpha, 
            v = 1))
    }
    else character()
}



# Read in the data
data<-read.csv("Data_post_review.csv")
head(data)
Treatments<-unique(data$treatment)
Treatments<-sort(Treatments)

# Checking a few distributions

par(mfrow=c(1,2))

# Note post-review valentina wants these surfaces on the raw scale
# I will fit on the log scale then exponentiate the predictions


hist(data$triglycerides)
hist(log(data$triglycerides))
# ln transform trig data - very skewed
data$triglycerides<-log(data$triglycerides)

hist(data$ucp1)
hist(log(data$ucp1))
# ln transform data - very skewed
data$ucp1<-log(data$ucp1)

hist(data$ucp1.batw)
hist(log(data$ucp1.batw))
# ln transform data - very skewed
data$ucp1.batw<-log(data$ucp1.batw)

hist(data$P4)
hist(log(data$P4))
# ln transform data - very skewed
data$p4<-log(data$P4)

hist(data$FSH)
hist(log(data$FSH))
# ln transform data - very skewed
data$FSH<-log(data$FSH)

hist(data$LH)
hist(log(data$LH))
# ln transform data - very skewed
data$LH<-log(data$LH)

hist(data$body.fat)
hist(log(data$body.fat))
# ln transform BF data - very skewed
data$body.fat<-log(data$body.fat)

hist(data$red.liver)
hist(log(data$red.liver))
data$red.liver<-log(data$red.liver)

# hmmm fgf data seem somehow zero infated
par(mfrow=c(2,3))
hist(data$fgf21)
hist(data$fgf21[which(data$treatment == "blank")])
hist(data$fgf21[which(data$treatment == "dht")])

# Could be that is multi-model for different diets - logging could be OK, so let's do it
hist(log(data$fgf21))
hist(log(data$fgf21[which(data$treatment == "blank")]))
hist(log(data$fgf21[which(data$treatment == "dht")]))

data$fgf21<-log(data$fgf21)

####################

# Note post reveiw we now need these surfaces

# Uterus weight in mg (was in gs)
hist(data$uterus.weight.mg)
# Seems OK

# GTT
hist(data$GTT_AUC)
# Maybe a bit skewed - or could be one outlier

################################################################
################### 3D Plots of the two groups #################
################################################################

pdf("3D.by.group_Rev.pdf", width=15, height=5)

par(mfrow=c(1,3), mar=c(5,5,5,5))

# How many levels and colors should there be on the surface
nlev<-8
no.cols<-256

# I will also create a list of the x and y variables in each set of predicted data
# I have given these very technical names as they will be aded to the plot, which
# theoretically should be publication quality
labs<-list()
labs[[1]]<-c("Protein eaten\n(kJ/mse/cage/d)", "Carbohydrate eaten\n(kJ/mse/cage/d)")
labs[[2]]<-c("Protein eaten\n(kJ/mse/cage/d)", "Fat eaten\n(kJ/mse/cage/d)")
labs[[3]]<-c("Carbohydrate eaten\n(kJ/mse/cage/d)", "Fat eaten\n(kJ/mse/cage/d)")

# This specifies the color scheme for surface
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# For the medians, note we use the median of all animals in the dataset.
PV<-list()

# Need to specify common limits on the polygons, and will also be used for the plots - range(data$intake.P), range(data$intake.C), range(data$intake.F)
p.limits<-c(0, 25)
c.limits<-c(0, 60)
f.limits<-c(0, 50)

fit.resolution<-101

# find the polygons for each treatment
for(i in 1:length(Treatments)){
	
	# Subset the data
	data.sub<-data[which(data$treatment == Treatments[i]),]
	
	# Get the polygons for each treatment
	Pred.Values<-list()
	Pred.Values[[1]]<-findConvex(x=data.sub$intake.P, y=data.sub$intake.C, rgnames=c("intake.P","intake.C"), res=fit.resolution, x.limits=p.limits, y.limits=c.limits)
	Pred.Values[[1]]$intake.F<-median(data$intake.F)
	Pred.Values[[2]]<-findConvex(x=data.sub$intake.P, y=data.sub$intake.F, rgnames=c("intake.P","intake.F"), res=fit.resolution, x.limits=p.limits, y.limits=f.limits)
	Pred.Values[[2]]$intake.C<-median(data$intake.C)
	Pred.Values[[3]]<-findConvex(x=data.sub$intake.C, y=data.sub$intake.F, rgnames=c("intake.C","intake.F"), res=fit.resolution, x.limits=c.limits, y.limits=f.limits)
	Pred.Values[[3]]$intake.P<-median(data$intake.P)
	
	# Save the polygons
	PV[[i]]<-Pred.Values

}

# Specify the model
model.use<-"outcome ~ s(intake.P, k=k1) + s(intake.C, k=k1) + s(intake.F, k=k1) + s(intake.P, intake.C, k=k2) + s(intake.P, intake.F, k=k2) + s(intake.C, intake.F, k=k2)"

write.table("3D GAM Results", file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F)

# Set upper DFs
k1<-4
k2<-6

# Outcomes That we are interested in for the paper
outcomes<-c("no.CL", "body.weight", "adiponectin", "cholesterol", "basal.glucose", "triglycerides", "GTT_AUC", "P4", "FSH", "LH")
families<-c(rep("gaussian", length(outcomes)))
families[1]<-"nb"

# Which of these have we log transformed
ln_transf<-rep(0, length(outcomes))
ln_transf[c(6, 8:10)]<-1

# Other possible outcomes of interest are: "body.fat", "uterus.weight", "adipo.size", "red.liver", "ucp1", "ucp1.batw", "fgf21", "no.cycles"


for(o in 1:length(outcomes)){
	
	# Specify the outcome
	data$outcome<-data[,outcomes[o]]
	
	# List to hold the surfaces
	fits<-list()
	mins<-NA
	maxs<-NA
	
	# Fit the GAMs and get predictions for the two treatment groups
	for(j in 1:length(Treatments)){
				
		# Subset the data and remove NAs
		data.sub<-data[which(data$treatment == Treatments[j]),]
		
		# Values for predictions
		Pred.Values<-PV[[j]]
				
		# Fit the GAM
		GAM<-gam(as.formula(model.use), data=data.sub, method="REML", family=families[o])
		res<-round(summary(GAM)$s.table, 3)
		res<-as.data.frame(cbind(row.names(res), res))
		names(res)[1]<-"Coef"
		
		write.table("", file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F, append=T)
		write.table(outcomes[o], file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F, append=T)
		write.table("", file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F, append=T)
		write.table(paste("Treatment = ", Treatments[j], sep=""), file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F, append=T)
		write.table(res, file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=names(res), append=T)
		
		misc<-paste0("n = ", summary(GAM)$n, ": % Dev Exp = ", summary(GAM)$dev.expl * 100)
		write.table(misc, file="3D_GAM_Results_Rev.csv", sep=",", row.names=F, col.names=F, append=T)
		
		# Save the models and datasets
		Fitted.Values<-list()
		for(k in 1:3){
			Fitted.Values[[k]]<-predict(GAM, newdata=Pred.Values[[k]], se.fit=T, type="response")			
		}
			
		fits[[j]]<-Fitted.Values
	
	}
		
	# Find the minimum and max fitted values
	vals<-c(unlist(fits[[1]][[1]]$fit), unlist(fits[[1]][[2]]$fit), unlist(fits[[1]][[3]]$fit), unlist(fits[[2]][[1]]$fit), unlist(fits[[2]][[2]]$fit), unlist(fits[[2]][[3]]$fit))
	
	if(ln_transf[o] == 1){
		vals<-exp(vals)
	}
	
	mn<-min(vals, na.rm=T)
	mx<-max(vals, na.rm=T)
	
	# Fit the GAMs and get predictions for the two treatment groups
	for(j in 1:length(Treatments)){
		
		# Values for predictions
		Pred.Values<-PV[[j]]
		Fitted.Values<- fits[[j]]
		
		for(i in 1:3){
		
			# Pretty comes up with nice values of x and y over which to fit
			if(i < 3){px<-p.limits}
			if(i == 3){px<-c.limits}
			if(i == 1){py<-c.limits}
			if(i > 1){py<-f.limits}
			
			# Uses px and py to generate the x and y axes
			x.new<-seq(min(px), max(px), len=fit.resolution)
			y.new<-seq(min(py), max(py), len=fit.resolution)
			
			surf<-matrix(Fitted.Values[[i]]$fit, nrow=sqrt(dim(Pred.Values[[i]])[1]))
			
			# Exponentiating the surface for the last 3 traits
			if(ln_transf[o] == 1){
				surf<-exp(surf)
			}
			
			locs<-round((range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols, 3)
		
			image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab=labs[[i]][1], ylab=labs[[i]][2], axes=FALSE)
			axis(1)
			axis(2)
			contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn, mx), nlev), labcex=1)
			if(i == 2){mtext(paste(Treatments[j], " - ", outcomes[o], sep=""), line=1, cex=1.5)}
	
		}
				
	}
	
	rm(fits)
	rm(vals)
	rm(mn)
	rm(mx)

}

dev.off()


rm(PV)

################################################################
####### 2D Plots of the PF surfaces with subtraction ###########
################################################################

pdf("3D.by.group.diff_Rev.pdf", width=15, height=5)

par(mfrow=c(1,3), mar=c(5,5,5,5))

# I will also create a list of the x and y variables in each set of predicted data
# I have given these very technical names as they will be aded to the plot, which
# theoretically should be publication quality
labs<-c("Carbohydrate eaten\n(kJ/mse/cage/d)", "Fat eaten\n(kJ/mse/cage/d)")


PV<-list()
# find the polygons for each treatment
for(i in 1:length(Treatments)){
	
	# Subset the data
	data.sub<-data[which(data$treatment == Treatments[i]),]
	
	# Get the polygons for each treatment
	Pred.Values<-list()
	Pred.Values[[1]]<-findConvex(x=data.sub$intake.C, y=data.sub$intake.F, rgnames=c("intake.C","intake.F"), res=fit.resolution, x.limits=c.limits, y.limits=f.limits)
	Pred.Values[[1]]$intake.P<-median(data$intake.P)
	
	# Save the polygons
	PV[[i]]<-Pred.Values

}

# List to hold the surfaces
fits<-list()
fits_se<-list()
mins<-NA
maxs<-NA

for(o in 1:length(outcomes)){

	# find the right outcome
	data$outcome<-data[,outcomes[o]]

	# Fit the GAMs and get predictions for the two treatment groups
	for(j in 1:length(Treatments)){
				
		# Subset the data and remove NAs
		data.sub<-data[which(data$treatment == Treatments[j]),]
		
		# Values for predictions
		Pred.Values<-PV[[j]][[1]]
				
		# Fit the GAM
		GAM<-gam(as.formula(model.use), data=data.sub, method="REML", family=families[o])
		res<-round(summary(GAM)$s.table, 3)
		res<-as.data.frame(cbind(row.names(res), res))
		names(res)[1]<-"Coef"
				
		# Save the models and datasets
		fits[[j]]<-predict(GAM, newdata=Pred.Values, se.fit=T, type="response")
		fits_se[[j]]<-predict(GAM, newdata=Pred.Values, se.fit=T, type="link")
		
	}
	
	# Find the minimum and max fitted values
	vals<-c(unlist(fits[[1]]$fit), unlist(fits[[2]]$fit))
	
	# Exponentiating for the last 3 traits
	if(ln_transf[o] == 1){
		vals<-exp(vals)
	}
	mn<-min(vals, na.rm=T)
	mx<-max(vals, na.rm=T)
	
	# Fit the GAMs and get predictions for the two treatment groups
	for(j in 1:length(Treatments)){
		
		# Values for predictions
		Pred.Values<-PV[[j]][[1]]
		Fitted.Values<-fits[[j]]
		# Exponentiating for the last 3 traits
		if(ln_transf[o] == 1){
			Fitted.Values$fit<-exp(Fitted.Values$fit)
		}
		
		px<-c.limits
		py<-f.limits
			
		# Uses px and py to generate the x and y axes
		x.new<-seq(min(px), max(px), len=fit.resolution)
		y.new<-seq(min(py), max(py), len=fit.resolution)
			
		surf<-matrix(Fitted.Values$fit, nrow=sqrt(dim(Pred.Values)[1]))
		locs<-round((range(unlist(Fitted.Values$fit), na.rm=TRUE) - mn) / (mx-mn) * no.cols, 3)
		image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab=labs[1], ylab=labs[2], axes=FALSE)
		axis(1)
		axis(2)
		contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn, mx), nlev), labcex=1)
		mtext(paste(Treatments[j], " - ", outcomes[o], sep=""), line=1, cex=1.5)
	
	}
				
	# Calculate the difference and SE of difference
	differences<-fits_se[[1]]$fit - fits_se[[2]]$fit
	se.differences<-sqrt(fits_se[[1]]$se.fit^2 + fits_se[[2]]$se.fit^2)
	signif.differences<-rep(NA, length(se.differences))
	
	# work out whether the difference has a CI spanning 0
	pos<-which(differences > 0)
	neg<-which(differences < 0)
	
	signif.differences[pos]<-as.numeric((differences[pos] - se.differences[pos] * 3.29053) >= 0)
	signif.differences[neg]<-as.numeric((differences[neg] + se.differences[neg] * 3.29053) <= 0)
			
	signif.differences[pos]<-signif.differences[pos] + as.numeric((differences[pos] - se.differences[pos] * 2.57583) >= 0)
	signif.differences[neg]<-signif.differences[neg] + as.numeric((differences[neg] + se.differences[neg] * 2.57583) <= 0)
						
	signif.differences[pos]<-signif.differences[pos] + as.numeric((differences[pos] - se.differences[pos] * 1.96) >= 0)
	signif.differences[neg]<-signif.differences[neg] + as.numeric((differences[neg] + se.differences[neg] * 1.96) <= 0)
	
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	mn<-0
	mx<-1
	
	px<-c.limits
	py<-f.limits
	
	# Colors to show significance
	signif.map<-c("lightgrey", "orchid1", "orchid3", "orchid4")
	
	# Uses px and py to generate the x and y axes
	x.new<-seq(min(px), max(px), len=fit.resolution)
	y.new<-seq(min(py), max(py), len=fit.resolution)
	
	surf<-matrix(signif.differences, nrow=sqrt(length(signif.differences)))
	up<-max(unlist(signif.differences), na.rm=TRUE)
	image(x.new, y.new, surf, col=signif.map[c(1:(up+1))], xlab=labs[1], ylab=labs[2], axes=FALSE)
	axis(1)
	axis(2)
	mtext("Significance of difference", line=1, cex=1.5)
	
	legend(37.5, 47.5, c("Non-Significant", "Significant at 0.05", "Significant at 0.01", "Significant at 0.001"), pch=16, col=signif.map, cex=1.1)
	

}

dev.off()


################################################################
################### MIXTURE MODELS FOR INTAKE ##################
################################################################

# File to write results to
res.file<-"Results_Mouse_PCOS_Intake_Rev.csv"

# Make sure proportions are closed to 1
sum<-data$diet.P + data$diet.C + data$diet.F
data$p_P<-data$diet.P/sum
data$p_C<-data$diet.C/sum
data$p_F<-data$diet.F/sum

# Open the pdf file for plotting
pdf("Intakes_Rev.pdf", height=5, width=10)

# Set the layout
par(mfrow=c(1,2), mar=c(5,5,5,1))	

# Set the resolution of the surface
surface.resolution<-501

# colour pallette
pall<-fields.colors(30)

# Color to show point for self-selected diet
point.col<-"black"

# How many values to round surface
round.surf<-3

# This specifies the color scheme for surface - it is actually a function that returns a function
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# How many different colours should we use on the plot
no.cols<-256

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# How many levels should there be on the surface
nlev<-5

treatments<-c("dht", "blank")

# Titles for each plot
titles<-list()
titles[[1]]<-"Blank"
titles[[2]]<-"DHT"

# Labels for each
labels<-c("Protein (%)", "Carbohydrate (%)", "Fat (%)")

# Find minimal and maximal values so as to scale sensibly
mn<-2
mx<-3.75

# Find minimal and maximal values so as to scale sensibly
#mn<-2
#mx<-4


## estimate convex hull and predict
mdff2<-findConvex.prop(data$p_P, data$p_C, c("p_P","p_C"), surface.resolution)
mdff2$p_F<-with(mdff2, 1 - p_P - p_C)

## plot the RMT surface
iso.lines<-seq(1, 0, -0.2)

for(k in 1:length(treatments)){
	
	# Find the subset for the treatment
	summary.stats<-data[which(data$treatment == treatments[k]),]
	
	# Create table for results
	if(k == 1){write.table("", file=res.file, sep=",", row.names=F, col.names=F)}
		
	# Variable to hold the four models of he scheffe's polynomials
	mmods<-list()
	
	# Fit the intercept model
	mmods[[1]]<-lm(ave.intake ~ 1, data=summary.stats)
	
	# Fit Scheffes polynomials over variance effect size for ith trait weighted by inverse sample variance - note not enough data for model 3
	for(i in 1:3){
		model<-i
		# Skip model 3
		if(i == 3){model<-4}
		# Fit the model
		mmods[[i+1]]<-MixModel(frame=summary.stats, response="ave.intake", mixcomps=c("p_P","p_C","p_F"), model=model)
	}
	
	# Find minimal model based on AIC
	AICs<-unlist(lapply(mmods, AIC))
	deltas<-AICs - min(AICs)
	options<-which(deltas <= 2)
	min.model<-min(options)
	model.AIC<-mmods[[min.model]]
	
	# Write the results to the table
	write.table("", file=res.file, sep=",", row.names=F, col.names=F, append=T)
	write.table(paste("Model ", min.model, " favoured by AIC."), file=res.file, sep=",", row.names=F, col.names=F, append=T)
	write.table(treatments[k], file=res.file, sep=",", row.names=F, col.names=F, append=T)
	res.k<-as.data.frame(round(summary(model.AIC)$coef[,c(1:3)], 4))
	res.k$df<-(dim(summary.stats)[1]) - (dim(res.k)[1])
	res.k$p<-(round(summary(model.AIC)$coef[,4], 4))
	res.k<-cbind(row.names(res.k), res.k)
	colnames(res.k)[1]<-"Coef."
	write.table(res.k, file=res.file, sep=",", row.names=F, col.names=colnames(res.k), append=T)
	
	# Get the predicted surface
	mdff2$fit<-predict(model.AIC, newdata=mdff2)
	surf<-matrix(mdff2$fit, nrow=sqrt(dim(mdff2)[1]))
	surf<-round(surf, round.surf)
	
	# Scale colors
	locs<-(range(surf, na.rm=TRUE) - mn) / (mx-mn) * no.cols
	
	# Actually plots the surface using all of the above info above
	plot(-10, -10, bty="n", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xlab="", ylab="", xaxt="n", yaxt="n")
	
	# Adds some axes
	axis(1, at = seq(0, 1, 0.2), labels=seq(0, 1, 0.2)*100)
	axis(2, at = seq(0, 1, 0.2), labels=seq(0, 1, 0.2)*100)
		
	# Add the Isolines
	for(i in 1:length(iso.lines)){
		abline(a = iso.lines[i], b=-1)
	}
	
	# Add the surface
	image(seq(0, 1, length.out = surface.resolution), seq(0, 1, length.out = surface.resolution), surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="", add=T)
	# Adds a contour over the top (can add a title using main)
	contour(seq(0, 1, length.out = surface.resolution), seq(0, 1, length.out = surface.resolution), surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
	
	# Add the axes labels
	mtext(labels[1], side=1, line=2.25)
	mtext(labels[2], side=2, line=2.25)
	text(0.55, 0.55, labels[3], srt=-45, cex=1.1)
	mtext(titles[[k]], cex=1.5, line=2)

}

dev.off()
