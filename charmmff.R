################################################################################
# The following is an exercise on reproducing step by step the elements which
# compose the CHARMM force-field. The energy expressions are very simple and
# could be written in a few lines of code, but the intention of this R
# workbook is to make it very context rich, including step by step and even
# redundant information.
#
# Written by: Mauricio Esguerra
# Date: November 26, 2013
################################################################################


################################################################################
# BONDED (INTRAMOLECULAR) INTERACTIONS
################################################################################
################################################################################
# BONDS AS HARMONIC POTENTIALS
################################################################################
#Force constants for carbon single, double, and triple bonds
k1=100; k2=200; k3=400

#Equilibrium distances for carbon single, double, and triple bonds in Ångstrom
b1=1.5; b2=1.3; b3=1.2

#Sampled variable
#b=runif(100, 0.5, 2.5)
b=seq(0.5, 2.5, ((2.5-0.5)/50))

#Harmonic Potential function
vbond1 <- k1*(b-b1)^2
vbond2 <- k2*(b-b2)^2
vbond3 <- k3*(b-b3)^2

################################################################################
# ANGLES AS HARMONIC POTENTIALS
################################################################################



################################################################################
# DIHEDRAL ANGLES AS FOURIER SERIES
################################################################################

k3=2.5; k2=5; k1=10
n3=3; n2=2; n1=1
delta=0
phi=seq(0, 2*pi, ((2*pi)/100))

vdihe1 <- k1*(1+(cos(n1*phi-delta)))
vdihe2 <- k2*(1+(cos(n2*phi-delta)))
vdihe3 <- k3*(1+(cos(n3*phi-delta)))
#vdihet <- vdihe1+vdihe2+vdihe3
#vdihe1 <- (sin(phi))*n1
#vdihe2 <- (sin(phi))*n2
#vdihe3 <- (sin(phi))*n3
vdihet <- vdihe1+vdihe2+vdihe3


################################################################################
# NON-BONDED (INTERMOLECULAR) INTERACTIONS
################################################################################












################################################################################
# POTENTIAL PLOTS
################################################################################
#pdf(file="charmm_potential.pdf",family="Helvetica", width=6,height=10)                           
par(mfcol=c(2,1))                                                                               
par(cex=0.4)                                                                                    
par(mar=c(4, 4, 0, 0), oma=c(6,6,7,2))                                                          
par(tcl=-0.25)                                                                                  
par(mgp=c(2,0.6,0)) 

xlim=range(0.5, 2.5) 
ylim=range(0,400)
plot(b, vbond1, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
lines(b, vbond2, col="red",type="o")
lines(b, vbond3, col="gray",type="o")
axis(1, cex.axis=1.0, at=seq(0, 2.5,0.5))
axis(2, cex.axis=1.0, at=seq(0, 400,50),las=1) 
box(which="plot",col="black")
legend("topright",col=c("blue","red","gray"),lty=1,legend=c("single","double","triple"), cex=0.8)


xlim=range(0, 360) 
ylim=range(0,30)
plot(phi*60, vdihe1, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
lines(phi*60, vdihe2, col="red",type="o")
lines(phi*60, vdihe3, col="gray",type="o")
lines(phi*60, vdihet, col="green",type="o")
axis(1, cex.axis=1.0, at=seq(0, 360,30))
axis(2, cex.axis=1.0, at=seq(0, 30,5),las=1) 
box(which="plot",col="black")
legend("topright",col=c("blue","red","gray","green"),lty=1,legend=c("n1","n2","n3","sum"), cex=0.8)

