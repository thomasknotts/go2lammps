__author__ = 'Anthony Gillespie'

import re
from GOLammpsInfo import GOLammpsInfo

class GOLammpsReader:

    # Utility variables
    currFileName = ""
    residuePat = "G\d+\s*"
    floatPat = "[0-9.-]+\s*"

    def __init__(self, info):
        self.info = info

    def read(self):

        self.read_xyz(self.info.go_pdb)
        self.read_num_chains(self.info.go_pdb)
        self.read_atom_masses(self.info.top)

        self.read_pair_coeffs(self.info.param)
        self.read_native_contacts(self.info.param)
        self.read_bond_coeffs(self.info.param)
        self.read_angle_coeffs(self.info.param)
        self.read_dihedral_coeffs(self.info.param)
        self.read_improper_coeffs(self.info.param)

        if self.info.original_pdb != "":
            self.read_amino_acid_identities(self.info.original_pdb)

        self.info.report_status()
        self.info.convert()


    def filter_lines(self, parsePattern, beginPattern, endPattern):
        lines = []
        file = open(self.currFileName)

        if beginPattern == '':
            foundBeginPattern = True
        else:
            foundBeginPattern = False

        for line in file:
                if foundBeginPattern and re.match(parsePattern, line):
                    lines.append(line)
                elif not foundBeginPattern and re.match(beginPattern, line):
                    foundBeginPattern = True
                elif endPattern != '' and re.match(endPattern, line):
                    break

        file.close()
        return lines

    def read_xyz(self, fileName):
        self.currFileName = fileName
        lines = self.filter_lines('.*ATOM.*','','')
        for line in lines:
            x = float(line.split()[5])
            y = float(line.split()[6])
            z = float(line.split()[7])
            self.info.coordinates.append((x,y,z))

    def read_num_chains(self, fileName):
        self.currFileName = fileName
        lines = self.filter_lines('.*ATOM.*','','')
        prev = 0
        totalParticleNum = 0
        for line in lines:
            chainParticleNum = int(line.split()[1])
            if chainParticleNum < prev:
                self.info.chainEnds.append(totalParticleNum)
            totalParticleNum += 1
            prev = chainParticleNum
        self.info.chainEnds.append(totalParticleNum)

    def read_atom_masses(self,fileName):
         self.currFileName = fileName
         lines = self.filter_lines('MASS','','')
         for line in lines:
             mass = int(round(float(line.split()[3])))
             if mass == 0:  # For some reason, all the histidines have a mass of 0
                 mass = 138
             self.info.masses.append(mass)

    def read_pair_coeffs(self,fileName):
        self.currFileName = fileName
        linePattern = self.residuePat + self.floatPat + self.floatPat
        lines = self.filter_lines(linePattern,'\s*NONBONDED.*','\s*NBFIX\s*')
        for line in lines:
            epsilon = float(line.split()[2])
            sigma = float(line.split()[3])
            self.info.nonbonded.append((epsilon, sigma))

    def read_native_contacts(self,fileName):
        self.currFileName = fileName
        linePattern = self.residuePat + self.residuePat + self.floatPat + self.floatPat
        lines = self.filter_lines(linePattern,'\s*NBFIX\s*','')
        for line in lines:
            residue1 = int((line.split()[0]).split('G')[1])
            residue2 = int((line.split()[1]).split('G')[1])
            epsilon = -1.0*float(line.split()[2])
            sigma = float(line.split()[3])
            self.info.nbfix.append((residue1, residue2, epsilon, sigma))

    def read_bond_coeffs(self,fileName):
        self.currFileName = fileName
        linePattern = self.residuePat + self.residuePat + self.floatPat + self.floatPat
        lines = self.filter_lines(linePattern,'\s*BOND\s*','\s*ANGLE\s*')
        for line in lines:
            residue1 = int((line.split()[0]).split('G')[1])
            residue2 = int((line.split()[1]).split('G')[1])
            K = float(line.split()[2])
            r0 = float(line.split()[3])
            self.info.bondCoeffs.append((residue1, residue2, K, r0))

    def read_angle_coeffs(self,fileName):
        self.currFileName = fileName
        linePattern = self.residuePat + self.residuePat + self.residuePat + self.floatPat + self.floatPat
        lines = self.filter_lines(linePattern,'\s*ANGLE\s*','\s*DIHEDRAL\s*')
        for line in lines:
            residue1 = int((line.split()[0]).split('G')[1])
            residue2 = int((line.split()[1]).split('G')[1])
            residue3 = int((line.split()[2]).split('G')[1])
            K = float(line.split()[3])
            theta0 = float(line.split()[4])
            self.info.angleCoeffs.append((residue1, residue2, residue3, K, theta0))

    def read_dihedral_coeffs(self,fileName):
        self.currFileName = fileName
        linePattern = self.residuePat + self.residuePat + self.residuePat + self.residuePat + self.floatPat + "[1-4]\s+" + self.floatPat
        lines = self.filter_lines(linePattern,'\s*DIHEDRAL\s*','\s*NONBONDED\s*')
        for line in lines:
            residue1 = int((line.split()[0]).split('G')[1])
            residue2 = int((line.split()[1]).split('G')[1])
            residue3 = int((line.split()[2]).split('G')[1])
            residue4 = int((line.split()[3]).split('G')[1])
            K = float(line.split()[4])
            n = int(line.split()[5])
            d = float(line.split()[6]) #Angle must be an int in LAMMPS
            weightingFactor = 0
            self.info.dihedralCoeffs.append((residue1, residue2, residue3, residue4, K, n, d, weightingFactor))

    def read_improper_coeffs(self,fileName):
        self.currFileName = fileName
        # The GO model builder does not include impropers

    def read_amino_acid_identities(self,fileName):
        self.currFileName = fileName
        #These are the hydropathy indices of the 20 naturally-ocurring amino acids
        aminoacids = dict( ALA = 1.8, ARG = -4.5, ASN = -3.5, ASP = -3.5, CYS = 2.5,
                            GLU = -3.5, GLN = -3.5, GLY = -0.4, HIS = -3.2, ILE = 4.5,
                            LEU = 3.8, LYS = -3.9, MET = 1.9, PHE = 2.8, PRO = 1.6,
                            SER = -0.8, THR = -0.7, TRP = -0.9, TYR = -1.3, VAL = 4.2 )
        linePattern = "\s*SEQRES\s*"
        lines = self.filter_lines(linePattern,'','')
        for line in lines:
            for string in line.split():
                if string in aminoacids:
                     self.info.aminoAcids.append(string)
