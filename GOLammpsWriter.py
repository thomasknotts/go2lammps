__author__ = 'Anthony Gillespie'

import time
from GOLammpsInfo import GOLammpsInfo

class GOLammpsWriter:

    def __init__(self, info):
        self.info = info

    def write(self):
        self.write_in_file()
        self.write_mod_file()
        self.write_data_file()

    def write_in_file(self):
        file = open(self.info.outputName + ".in",'w')
        file.write("""# Created on: %s

units             real
neigh_modify      delay 2 every 1

atom_style        molecular
bond_style        harmonic
angle_style       charmm
dihedral_style    eten
pair_style        lj/cut/eten 395.5
pair_modify       mix arithmetic
pair_modify       shift yes         # Turns on shifting in lj energies
special_bonds lj  0.0 0.0 1.0       # Calculate forces between particles in 1-4 position
neighbor          20.0 bin

# READ IN DATA
read_data         %s.data
include           %s_native_contacts.mod

# SET UP simulation parameters
variable          t world 300.0
variable          n world 1
fix               myfix all nvt temp $t $t 100.0 tchain 3
velocity          all create 300 12345678 dist uniform
timestep          3

# OUTPUTS
thermo            1000
thermo_style      custom step ebond eangle edihed eimp evdwl pe ke etotal
dump              1 all xyz 1000 %s.xyz

# RUN
run               50000

        """ % (time.strftime("%c"),self.info.outputName,self.info.outputName,self.info.outputName))
        file.close()


    def write_mod_file(self):
        file = open(self.info.outputName + "_native_contacts.mod",'w')
        for pair_coeff in self.info.nbfix:
            residue1 = pair_coeff[0]
            residue2 = pair_coeff[1]
            epsilon = pair_coeff[2]
            sigma = pair_coeff[3]
            line = "pair_coeff     %d      %d     %f     %f\n" % (residue1, residue2, epsilon, sigma)
            file.write(line)
        file.close()


    def write_data_file(self):
        file = open(self.info.outputName + ".data",'w')
        self.write_header(file)
        self.write_masses(file)
        self.write_pair_coeffs(file)
        self.write_atom_coordinates(file)
        self.write_bond_coeffs(file)
        self.write_bonds(file)
        self.write_angle_coeffs(file)
        self.write_angles(file)
        self.write_dihedral_coeffs(file)
        self.write_dihedrals(file)
        file.close()

    def write_header(self,file):
        numAtoms = len(self.info.coordinates)
        numBonds = len(self.info.bondCoeffs)
        numAngles = len(self.info.angleCoeffs)
        numDihedrals = len(self.info.dihedralCoeffs)
        numImpropers = 0
        maxVal = self.info.get_max_abs_xyz()
        max = round(maxVal*5)
        if max < 791:
          max = 791
        min = max*-1
        file.write("""Created by go2lammps.py on: %s

        %d          atoms
        %d          bonds
        %d          angles
        %d          dihedrals
        %d          impropers

        %d          atom types
        %d          bond types
        %d          angle types
        %d          dihedral types
        %d          improper types

        %f      %f      xlo xhi
        %f      %f      ylo yhi
        %f      %f      zlo zhi
        """ % (time.strftime("%c"),
                numAtoms,numBonds,numAngles,numDihedrals,numImpropers,
                numAtoms,numBonds,numAngles,numDihedrals,numImpropers,
                min,max,min,max,min,max))

    def write_masses(self,file):
        file.write("\nMasses\n\n")
        count = 0
        for mass in self.info.masses:
            count += 1
            residue = str(count).ljust(5)
            mass = str(mass).ljust(5)
            line = "    %s     %s\n" % (residue, mass)
            file.write(line)

    def write_pair_coeffs(self,file):
        file.write("\nPair Coeffs\n\n")
        count = 0
        for nonbondedTuple in self.info.nonbonded:
            count += 1
            residueNum = str(count).ljust(5)
            epsilon = str(nonbondedTuple[0]).ljust(5)
            sigma = str(nonbondedTuple[1]).ljust(5)
            line = "    %s    %s    %s\n" % (residueNum, epsilon, sigma)
            file.write(line)

    def write_atom_coordinates(self,file):
        file.write("\nAtoms\n\n")
        for i,(x,y,z) in enumerate(self.info.coordinates):
            atomNum = str(i+1).ljust(5)
            moleculeNum = str(self.info.get_molecule_number(i+1)).ljust(5)
            atomType = atomNum
            x = str(x).ljust(7)
            y = str(y).ljust(7)
            z = str(z).ljust(7)
            line = "    %s    %s    %s    %s    %s    %s\n" % (atomNum, moleculeNum, atomType, x, y, z)
            file.write(line)

    def write_bond_coeffs(self,file):
        file.write("\nBond Coeffs\n\n")
        for i,(residue1, residue2, K, r0) in enumerate(self.info.bondCoeffs):
            bondType = str(i+1).ljust(5)
            K = str(int(K)).ljust(5) # K has to be an int in lammps
            r0 = str(r0).ljust(5)
            line = "    %s    %s    %s\n"  % (bondType,K,r0)
            file.write(line)

    def write_bonds(self,file):
        file.write("\nBonds\n\n")
        for i,(residue1, residue2, K, r0) in enumerate(self.info.bondCoeffs):
            identifier = str(i+1).ljust(5)
            bondType = identifier
            residue1 = str(residue1).ljust(5)
            residue2 = str(residue2).ljust(5)
            line = "    %s    %s    %s    %s\n" % (identifier,bondType,residue1,residue2)
            file.write(line)

    def write_angle_coeffs(self,file):
        file.write("\nAngle Coeffs\n\n")
        for i,(residue1, residue2, residue3, K, theta0) in enumerate(self.info.angleCoeffs):
            angleType = str(i+1).ljust(5)
            K = str(K).ljust(5)
            theta0 = str(theta0).ljust(12)
            line = "    %s    %s    %s    0    0\n" % (angleType, K , theta0)
            file.write(line)

    def write_angles(self,file):
        file.write("\nAngles\n\n")
        for i,(residue1, residue2, residue3, K, theta0) in enumerate(self.info.angleCoeffs):
            identifier = str(i+1).ljust(5)
            angleType = identifier
            residue1 = str(residue1).ljust(5)
            residue2 = str(residue2).ljust(5)
            residue3 = str(residue3).ljust(5)
            line = "    %s    %s    %s    %s    %s\n" % (identifier,angleType,residue1,residue2,residue3)
            file.write(line)

    def write_dihedral_coeffs(self,file):
        file.write("\nDihedral Coeffs\n\n")
        for i,(residue1, residue2, residue3, residue4, K, n, d, weightingFactor) in enumerate(self.info.dihedralCoeffs):
            dihedralType = str(i+1).ljust(5)
            K = str(K).ljust(12)
            n = str(n).ljust(3)
            d = str(d).ljust(5)
            weightingFactor = str(weightingFactor).ljust(5)
            line = "    %s    %s    %s    %s    %s\n" % (dihedralType, K, n, d, weightingFactor)
            file.write(line)

    def write_dihedrals(self,file):
        file.write("\nDihedrals\n\n")
        for i,(residue1, residue2, residue3, residue4, K, n, d, weightingFactor) in enumerate(self.info.dihedralCoeffs):
            identifier = str(i+1).ljust(5)
            dihedralType = identifier
            residue1 = str(residue1).ljust(5)
            residue2 = str(residue2).ljust(5)
            residue3 = str(residue3).ljust(5)
            residue4 = str(residue4).ljust(5)
            line = "    %s    %s    %s    %s    %s    %s\n" % (identifier,dihedralType,residue1,residue2,residue3,residue4)
            file.write(line)
