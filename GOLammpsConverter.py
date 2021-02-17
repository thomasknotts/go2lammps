__author__ = 'Anthony Gillespie'

from GOLammpsInfo import GOLammpsInfo
from GOLammpsReader import GOLammpsReader
from GOLammpsWriter import GOLammpsWriter

#Options
GOModel = True

if GOModel:

    # REPLACE THE FOLLOWING VARIABLES WITH 
    # THE NAMES OF THE FILES YOU ARE USING

    go_pdb = "GO_1BDD10-55.pdb"      # The .pdb file that came from the MMTSB server
    top = "GO_1BDD10-55.top"         # The .top file that came from the MMTSB server
    param = "GO_1BDD10-55.param"     # The .param file that came from the MMTSB server
    original_pdb = "1BDD10-55.pdb"       # The original PDB file you submitted to the MMTSB server
    outputName = "1bdd"         # The name you want the LAMMPS files to have

    info = GOLammpsInfo(go_pdb,top,param,original_pdb,outputName)
    reader = GOLammpsReader(info)
    writer = GOLammpsWriter(info)

    reader.read()
    writer.write()
