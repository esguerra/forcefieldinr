################################################################################
# Written by: Mauricio Esguerra
# Date: November 26, 2013
# Stockholm and Uppsala, Sweden
#
# The following is an exercise on reproducing step by step the elements which
# compose the CHARMM force-field. The energy expressions are very simple and
# could be written in a few lines of code, but the intention of this R
# workbook is to make it very context rich, including step by step and even
# redundant information.
#
# In CHARMM there are various parameter files containing the force constants and
# equilibrium values for bonds, angles, and dihedrals. For this example we are
# going to use the values found in the parameter file par_all36_prot.prm and 
# toppar_all36_prot_fluoro_alkanes.str. Relevant lines are extracted in the 
# following text.
#
# We have selected ethane as a test molecule.
# The cartoon corresponding to the topology of ethane is:
#
#    HA3        HA3
#      \       /
#   HA3-CT3---CT3-HA3
#      /       \
#    HA3        HA3
#
#Resi ETHA  0.0 ! ethane
#Group 
#Atom h13 ha3   0.09
#Atom c1  ct3  -0.27
#Atom c2  ct3  -0.27
#Atom h21 ha3   0.09
#Atom h22 ha3   0.09
#Atom h23 ha3   0.09
#Atom h11 ha3   0.09
#Atom h12 ha3   0.09
#Bond c1 h11  c1 h12   c1 h13
#Bond c1 c2
#Bond c2 h21  c2 h22   c2 h23
#ic h21 c2 c1  h13  0.00 0.00  180.0 0.0  0.0
#ic h22 c2 c1  h13  0.00 0.00   60.0 0.0  0.0
#ic h23 c2 c1  h13  0.00 0.00  -60.0 0.0  0.0
#ic h11 c1 c2  h21  0.00 0.00   60.0 0.0  0.0
#ic h12 c1 c2  h21  0.00 0.00  -60.0 0.0  0.0
#ic c1  c2 h21 h22  0.00 0.00    0.0 0.0  0.0
#patc firs none last none
#
#
# MASS    45 HA3    1.00800 ! alkane, CH3, new LJ params (see toppar_all22_prot_aliphatic_c27.str)
# MASS    52 CT3   12.01100 ! aliphatic sp3 C for CH3
# CT3  CT3   222.500     1.5300 ! ALLOW   ALI
#                ! alkane update, adm jr., 3/2/92
# HA3  CT3   322.000     1.1110 ! ALLOW   ALI
# ! alkane update, adm jr., 3/2/92
#
# DIHEDRALS
# !
# !V(dihedral) = Kchi(1 + cos(n(chi) - delta))
# !
# !Kchi: kcal/mole
# !n: multiplicity
# !delta: degrees
# !
# !atom types             Kchi    n   delta
# HA3  CT3  CT3   37.500    110.10   22.53   2.17900 ! ALLOW   ALI
# ! alkane update, adm jr., 3/2/92
# HA3  CT3  HA3   35.500    108.40    5.40   1.80200 ! ALLOW   ALI
# ! alkane update, adm jr., 3/2/92
# X    CT3  CT3  X        0.1525  3     0.00 ! ALLOW   ALI
# ! alkane, 4/98, yin and mackerell
# 
# !V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
# !
# !epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
# !Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
# !
# !atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
# !
# !carbons
# CT3      0.0       -0.0780    2.040 0.0 -0.01 1.9 ! alkane, 4/98, yin, adm jr.
# !hydrogens
# HA     0.000000  -0.022000     1.320000 ! ALLOW PEP ALI POL SUL ARO PRO ALC
# ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92
################################################################################





################################################################################
# BONDED (INTRAMOLECULAR) INTERACTIONS
################################################################################
################################################################################
# BONDS AS HARMONIC POTENTIALS
################################################################################
#Force constants for carbon single, double, and triple bonds
k1=222.50; k2=300; k3=450

#Equilibrium distances for carbon single, double, and triple bonds in Ångstrom
b1=1.5300; b2=1.3; b3=1.0

#Sampled variable
#b=runif(100, 0.5, 2.5)
b=seq(0.0, 2.5, ((2.5-0.0)/50))

#Harmonic Potential function
vbond1 <- k1*(b-b1)^2
vbond2 <- k2*(b-b2)^2
vbond3 <- k3*(b-b3)^2

################################################################################
# ANGLES AS HARMONIC POTENTIALS
################################################################################
#Force constants for carbon-carbon-carbon carbon-carbon-hydrogen
#k1=50; k2=100; k3=200

#Equilibrium angles
#b1=pi/2; b2=pi/3; b3=pi/4

#Sampled variable
#b=runif(100, 0.5, 2.5)
#b=seq(0.0, pi, ((pi-0.0)/50))

#Harmonic Potential function
#vangle1 <- k1*(b-b1)^2
#vangle2 <- k2*(b-b2)^2
#vangle3 <- k3*(b-b3)^2



################################################################################
# DIHEDRAL ANGLES AS COSINES OR FOURIER SERIES
################################################################################

k3=0.1525; k2=0.3; k1=0.6 # k=Vo/2, so the barrier for rotation is Vo.
n3=3; n2=2; n1=1    # number of barriers or minima.
delta=0
phi=seq(0, 2*pi, ((2*pi)/100))

#angles are given in radians so have to convert on ploting to degrees.
vdihe1 <- k1*(1+(cos(n1*phi-delta)))
vdihe2 <- k2*(1+(cos(n2*phi-delta)))
vdihe3 <- k3*(1+(cos(n3*phi-delta)))

vdihet <- vdihe1+vdihe2+vdihe3


################################################################################
# NON-BONDED (INTERMOLECULAR) INTERACTIONS
################################################################################
#lennard-jones potential
epsiloni=-0.0780 #kcal/mole * 10
epsilonj=-0.0220 #kcal/mole * 10
epsilonij=sqrt(epsiloni*epsilonj) #Lorentz-Berthelot
rmini=2.040*2
rminj=1.32*2
rminij=((rmini+rminj)/2)
rij=seq(0.0, 7.0, ((7.0-0.0)/100))
lj=epsilonij*( ((rminij/rij)**12) - (2*(rminij/rij)**6) )

#electrostatic potential





################################################################################
# POTENTIAL PLOTS
################################################################################
pdf(file="charmm_potential.pdf",family="Helvetica", width=6,height=10)                           
par(mfcol=c(3,1))                                                                               
par(cex=0.4)                                                                                    
par(mar=c(4, 4, 0, 0), oma=c(6,6,7,2))                                                          
par(tcl=-0.25)                                                                                  
par(mgp=c(2,0.6,0)) 

xlim=range(0.0, 2.5) 
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
ylim=range(0,1.4)
plot(phi*(360/(2*pi)), vdihe1, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
lines(phi*(360/(2*pi)), vdihe2, col="red",type="o")
lines(phi*(360/(2*pi)), vdihe3, col="gray",type="o")
#lines(phi*(360/(2*pi)), vdihet, col="green",type="o")
axis(1, cex.axis=1.0, at=seq(0, 360,30))
axis(2, cex.axis=1.0, at=seq(0, 2,0.2),las=1) 
box(which="plot",col="black")
legend("topright",col=c("blue","red","gray","green"),lty=1,legend=c("n1","n2","n3","sum"), cex=0.8)


xlim=range(2.5, 6.0) 
ylim=range(-0.1,0.1)
plot(rij,lj, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
#lines(b, vbond2, col="red",type="o")
#lines(b, vbond3, col="gray",type="o")
abline(v=rminij,col=3,lty=3)
axis(1, cex.axis=1.0, at=seq(2.5, 6.0, 0.5))
axis(2, cex.axis=1.0, at=seq(-0.1, 0.1, 0.05),las=1) 
box(which="plot",col="black")
#legend("topright",col=c("blue","red","gray"),lty=1,legend=c("single","double","triple"), cex=0.8)




dev.off()