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
# going to use the values found in the PARAMETER file par_all36_prot.prm and 
# TOPOLOGY file stream/prot/toppar_all36_prot_model.str. Relevant lines are extracted in the 
# following text.
#
#Furthermore in CHARMM a distinction is made between TOPOLOGY files, that is,
#those describing the connectivity of atoms in the molecule, and PARAMETER files
#those describing the parameters used in the classical potential expressions to
#describe bonds and angles.
#
# We have selected ethane as a test molecule.
#
#TOPOLOGY
#============================
# Resi ETHA        0.00 ! ethane, S. Fischer
# Group
# Group
# Atom h11 ha3     0.09 !  H11      H21
# Atom h12 ha3     0.09 !    \      /
# Atom h13 ha3     0.09 ! H12-C1--C2-H22
# Atom c1  ct3    -0.27 !    /      \
# Group                 !   H13       H23
# Atom h21 ha3     0.09
# Atom h22 ha3     0.09
# Atom h23 ha3     0.09
# Atom c2  ct3    -0.27
# Bond c1 h11  c1 h12  c1 h13
# Bond c1 c2   c2 h21  c2 h22  c2 h23
# ic  h11 c1 c2  h21  0.00 0.00    0.0 0.0  0.0
# ic  h11 c1 c2  h22  0.00 0.00  120.0 0.0  0.0
# ic  h11 c1 c2  h23  0.00 0.00  240.0 0.0  0.0
# ic  h12 c1 c2  h23  0.00 0.00  120.0 0.0  0.0
# ic  h13 c1 c2  h23  0.00 0.00  240.0 0.0  0.0

#PARAMETERS
#============================
# ATOMS
# MASS    45 HA3    1.00800 ! alkane, CH3, new LJ params (see toppar_all22_prot_aliphatic_c27.str)
# MASS    52 CT3   12.01100 ! aliphatic sp3 C for CH3
# 
# BONDS
# !
# !V(bond) = Kb(b - b0)**2
# !
# !Kb: kcal/mole/A**2
# !b0: A
# !
# !atom type Kb          b0
# CT3  CT3   222.500     1.5300 ! ALLOW   ALI
#                ! alkane update, adm jr., 3/2/92
# HA3  CT3   322.000     1.1110 ! ALLOW   ALI
# ! alkane update, adm jr., 3/2/92
#
#
# ANGLES
# !
# !V(angle) = Ktheta(Theta - Theta0)**2
# !
# !V(Urey-Bradley) = Kub(S - S0)**2
# !
# !Ktheta: kcal/mole/rad**2
# !Theta0: degrees
# !Kub: kcal/mole/A**2 (Urey-Bradley)
# !S0: A
# !
# !atom types     Ktheta    Theta0   Kub     S0
# HA3  CT3  CT3   37.500    110.10   22.53   2.17900 ! ALLOW   ALI
# ! alkane update, adm jr., 3/2/92
# HA3  CT3  HA3   35.500    108.40    5.40   1.80200 ! ALLOW   ALI


# DIHEDRALS
# !
# !V(dihedral) = Kchi(1 + cos(n(chi) - delta))
# !
# !Kchi: kcal/mole
# !n: multiplicity
# !delta: degrees
# !
# !atom types             Kchi    n   delta
# ! alkane update, adm jr., 3/2/92
# X    CT3  CT3  X        0.1525  3     0.00 ! ALLOW   ALI
# ! alkane, 4/98, yin and mackerell
# 

# IMPROPER
# !
# !V(improper) = Kpsi(psi - psi0)**2
# !
# !Kpsi: kcal/mole/rad**2
# !psi0: degrees
# !note that the second column of numbers (0) is ignored
# !
# !atom types           Kpsi                   psi0
# !

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

#Equilibrium distances for carbon single, double, and triple bonds in Ångström
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

#k3=0.1525; k2=0.3; k1=0.6 # k=Vo/2, so the barrier for rotation is Vo.
k3=0.1525; k2=k3; k1=k3
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
epsiloni=-0.0780 #kcal/mole
epsilonj=-0.0220 #kcal/mole
epsilonij=sqrt(epsiloni*epsilonj) #Lorentz-Berthelot Lorentz=geometric
rmini=2.040*2   #Ångström
rminj=1.32*2    #Ångström
rminij=((rmini+rminj)/2)          #Lorentz-Berthelot Berthelot=aritmetic
rij=seq(0.0, 7.0, ((7.0-0.0)/100))
Aij= ((rminij/rij)**12)
Bij= (2*(rminij/rij)**6)
lj=epsilonij*( Aij - Bij)

#electrostatic potential
# It is Coulomb's law in AKMA units, ie length in A, energy in kcal/mol and charge in |e|.
# A relative dielectric permittivity (dimensionless) is allowed for in the denominator 
# (keyword EPS, default value 1.0). Numerically this boils down to:
# E= 332*(qi*qj)/(rij*eps)
# Kcoul = 1/(4*pi*epsilon-zero)
# epsilon-zero = 8.8542 x 10^-12 C^2 m^-1 J^-1
# 1 m = 10^10 A
# 1 C = 1.6022 x 10 ^19 esu
# 1 J = 6.0221 x 10 23 / 4184 kcal mol^-1
# Kcoul = 1 J m / 4*pi*(8.8542 x 1o ^-12)C^2 == 332 kcal * mol^-1 * esu^-2
qi=-0.27  # carbon charge
qj=0.09   # hydrogen charge
e=1 # dielectric constant 1 for vacuum, 80 for water, 4 for protein interior 
rij=seq(0.0, 7.0, ((7.0-0.0)/100))
Kcoul=332
velecc= Kcoul*(qi*qi)/(e*rij)
velech= Kcoul*(qi*qj)/(e*rij)

################################################################################
# POTENTIAL PLOTS
################################################################################
pdf(file="charmm_potential.pdf",family="Helvetica", width=10,height=6)                           
par(mfcol=c(2,2))                                                                               
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
lines(phi*(360/(2*pi)), vdihet, col="green",type="o")
axis(1, cex.axis=1.0, at=seq(0, 360,30))
axis(2, cex.axis=1.0, at=seq(0, 2,0.2),las=1) 
box(which="plot",col="black")
legend("topright",col=c("blue","red","gray","green"),lty=1,legend=c("n1","n2","n3","sum"), cex=0.8)


xlim=range(2.5, 6.0) 
ylim=range(-0.2,0.2)
plot(rij,epsilonij*Aij, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
lines(rij,-epsilonij*Bij, pch=23, type="o", col="red", cex=0.8,
      xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
      panel.first = grid(),  xaxs="i", yaxs="i")
lines(rij,lj, pch=23, type="o", col="gray", cex=0.8,
      xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
      panel.first = grid(), xaxs="i", yaxs="i")
abline(v=rminij,col="gray",lty=2, lwd=2)
axis(1, cex.axis=1.0, at=seq(2.5, 6.0, 0.5))
axis(2, cex.axis=1.0, at=seq(-0.2, 0.2, 0.05),las=1) 
legend("topright",col=c("blue","red","gray"),lty=1,legend=c("repulsive","attractive","L-J"), cex=0.8)
box(which="plot",col="black")

xlim=range(0, 6.0) 
ylim=range(-20,20)
plot(rij,velecc, pch=23, type="o", col="blue", cex=0.8,
     xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
     panel.first = grid(), axes=FALSE, xaxs="i", yaxs="i")
lines(rij,velech, pch=23, type="o", col="red", cex=0.8,
      xlab = " ", main="", xlim = xlim, ylim = ylim, ylab= " ",
      panel.first = grid(),  xaxs="i", yaxs="i")
#abline(v=rminij,col="gray",lty=2, lwd=2)
axis(1, cex.axis=1.0, at=seq(0.5, 6.0, 0.5))
axis(2, cex.axis=1.0, at=seq(-20, 20, 5),las=1) 
legend("topright",col=c("blue","red"),lty=1,legend=c("C-C","C-H"), cex=0.8)
box(which="plot",col="black")


dev.off()