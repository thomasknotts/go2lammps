__author__ = 'Anthony Gillespie'

class GOLammpsInfo:

    go_pdb = ""
    top = ""
    param = ""
    original_pdb = ""
    outputName = ""
    converted = False

    coordinates = []
    chainEnds = []
    masses = []
    nonbonded = []
    nbfix = []
    bondCoeffs = []
    angleCoeffs = []
    dihedralCoeffs = []
    aminoAcids = []

    def __init__(self, go_pdb, top, param, original_pdb, outputName):
        self.go_pdb = go_pdb
        self.top = top
        self.param = param
        self.original_pdb = original_pdb
        self.outputName = outputName
        self.converted = False

    def get_min_max_xyz(self):
        length = len(self.coordinates)
        x_min = sorted(coord[0] for coord in self.coordinates)[0]
        x_max = sorted(coord[0] for coord in self.coordinates)[length-1]
        y_min = sorted(coord[1] for coord in self.coordinates)[0]
        y_max = sorted(coord[1] for coord in self.coordinates)[length-1]
        z_min = sorted(coord[2] for coord in self.coordinates)[0]
        z_max = sorted(coord[2] for coord in self.coordinates)[length-1]
        return (x_min, x_max, y_min, y_max, z_min, z_max)

    def get_max_abs_xyz(self):
        length = len(self.coordinates)
        maximums = []
        maximums.append(abs(sorted(coord[0] for coord in self.coordinates)[0]))
        maximums.append(abs(sorted(coord[0] for coord in self.coordinates)[length-1]))
        maximums.append(abs(sorted(coord[1] for coord in self.coordinates)[0]))
        maximums.append(abs(sorted(coord[1] for coord in self.coordinates)[length-1]))
        maximums.append(abs(sorted(coord[2] for coord in self.coordinates)[0]))
        maximums.append(abs(sorted(coord[2] for coord in self.coordinates)[length-1]))
        return sorted(maximums)[5]

    def report_status(self):
        print("len(coordinates): %d" % len(self.coordinates))
        print("len(chainEnds): %d" % len(self.chainEnds))
        print("len(masses): %d" % len(self.masses))
        print("len(nonbonded): %d" % len(self.nonbonded))
        print("len(nbfix): %d" % len(self.nbfix))
        print("len(bondCoeffs): %d" % len(self.bondCoeffs))
        print("len(angleCoeffs): %d" % len(self.angleCoeffs))
        print("len(dihedralCoeffs): %d" % len(self.dihedralCoeffs))
        print("len(aminoAcids): %d" % len(self.aminoAcids))

    def get_molecule_number(self,atomNum):
        for i,chainEnd in enumerate(self.chainEnds):
            if atomNum <= chainEnd:
                return i + 1

    def amino_acid_compare(self):

        aminoacidsMasses = dict( ALA = 71, ARG = 157, ASN = 114, ASP = 114, CYS = 103,
                            GLU = 128, GLN = 128, GLY = 57, HIS = 138, ILE = 113,
                            LEU = 113, LYS = 128, MET = 131, PHE = 147, PRO = 97,
                            SER = 87, THR = 101, TRP = 186, TYR = 163, VAL = 99 )

        massCount = 0
        for aminoAcid in self.aminoAcids:
            if self.masses[massCount] != aminoacidsMasses[aminoAcid]:
                #print("massCount: %d   mass %d   aminoacidMass %f   aminoAcid = %s" % (massCount, self.masses[massCount], aminoacidsMasses[aminoAcid], aminoAcid))
                print("%d  %d   %f   %s  ****************" % (massCount, self.masses[massCount], aminoacidsMasses[aminoAcid], aminoAcid))
            else:
                print("%d  %d   %f   %s" % (massCount, self.masses[massCount], aminoacidsMasses[aminoAcid], aminoAcid))
                massCount += 1

    def convert(self):
        self.convert_coordinates()
        self.convert_masses()
        self.convert_nonbonded()



    def convert_coordinates(self):
        return True
        #can shift coordinates to here if desired

    def convert_masses(self):
        for mass in self.masses:
            if mass == 0:   # For some reason, all the histidines have a mass of 0
                mass = 138

    def convert_nonbonded(self):
        convertedNonbonded = []
        for nonbondedTuple in self.nonbonded:
            epsilon = abs(nonbondedTuple[0]) #Sigma needs to be positive
            sigma = nonbondedTuple[1]*2 #Charmm stores epsilon/2 not espsilon
            convertedNonbonded.append((epsilon, sigma))
        self.nonbonded = convertedNonbonded

    def set_go_pdb(self, fileName):
        self.go_pdb = fileName

    def set_go_top(self, fileName):
        self.top = fileName

    def set_go_param(self, fileName):
        self.param = fileName

    def set_original_pdb(self, fileName):
        self.original_pdb = fileName

    def set_output_name(self, outputName):
        self.outputName = outputName
