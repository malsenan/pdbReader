from Bio.PDB import *
import math

#https://www.tutorialspoint.com/biopython/biopython_pdb_module.htm
#nodes for each atom in linked list w/ the atom name, coordinates, and next atom node
class AtomNode:

    def __init__(self, name = None, coords = None, residue = None):
        self.name = name    #string
        self.coordinates = coords  #list/array
        self.residue = residue
        self.nxt = None

    def getName(self):
        return self.name

    def getCoords(self):
        return self.coordinates

    def getRes(self):
        return self.residue

    def setName(self, newName):
        self.name = newName

    def setCoords(self, newCoords):
        self.coordinates = newCoords

    def setRes(self, newRes):
        self.residue = newRes

#list filled with atom nodes and statistics in a hashmap (dictionary)
class AtomList:

    def __init__(self):
        self.head = AtomNode()    #head of list
        self.stats = dict()     #hashmap with averages, standard deviations, etc

    #1
    def readPDB(self, protein):
        pdb = PDBList()
        pdb.retrieve_pdb_file(protein, pdir = '.', file_format = 'pdb')
        parser = PDBParser(PERMISSIVE = True, QUIET = True)
        data = parser.get_structure(protein, 'pdb' + protein + '.ent')
        models = list(data.get_models())
        chains = list(models[0].get_chains())
        allAtoms = list(chains[0].get_atoms())

        Nnum = 1
        CAnum = 1
        Cnum = 1
        curr = self.head
        for atom in allAtoms:
            if atom.get_name() == 'N':
                curr.setName(atom.get_name())
                curr.setCoords(atom.get_coord())
                curr.setRes(Nnum)
                Nnum += 1
                curr.nxt = AtomNode()
                curr = curr.nxt
            elif atom.get_name() == 'CA':
                curr.setName(atom.get_name())
                curr.setCoords(atom.get_coord())
                curr.setRes(CAnum)
                CAnum += 1
                curr.nxt = AtomNode()
                curr = curr.nxt
            elif atom.get_name() == 'C':
                curr.setName(atom.get_name())
                curr.setCoords(atom.get_coord())
                curr.setRes(Cnum)
                Cnum += 1
                curr.nxt = AtomNode()
                curr = curr.nxt

    #2
    def getStats(self):
        NtoCA = []  #array of all distances (for standard deviation)
        NtoCAtotal = 0  #total of distances between N and CA atoms

        CAtoC = []  #array of "   "
        CAtoCtotal = 0  #total of "   "

        CtoN = []   #array of "   "
        CtoNtotal = 0   #total of "   "

        angsNtoCAtoC = []     #array of angles of N to CA to C (in radians)
        angtotalNtoCAtoC = 0    #total sum of angles of N to CA to C

        angsCAtoCtoN = []
        angtotalCAtoCtoN = 0

        angsCtoNtoCA = []
        angtotalCtoNtoCA = 0

        v1 = []     #first vector that will be used to get angles ("previous vector")

        CAtoCA = []     #list of all distances between adjacent CA atoms
        CAtoCAtotal = 0     #sum of distances between adjacent CA atoms
        CA1 = []    #previous CA atom coordinate to find distance to current CA atom

        norm1 = []    #previous cross product of 2 vectors to get normal of plane to find torsional angles in residue 30

        curr = self.head
        while curr.getName() is not None:
            nextAtom = curr.nxt
            if nextAtom.getName() is not None:
                #distance formula for points in 3D space
                distance = math.sqrt((nextAtom.getCoords()[0] - curr.getCoords()[0]) ** 2 + (nextAtom.getCoords()[1] - curr.getCoords()[1]) ** 2 + (nextAtom.getCoords()[2] - curr.getCoords()[2]) ** 2)
                v2 = [nextAtom.getCoords()[0] - curr.getCoords()[0], nextAtom.getCoords()[1] - curr.getCoords()[1], nextAtom.getCoords()[2] - curr.getCoords()[2]]  # "current vector" to get bond angle
                if curr.getName() == 'N':
                    NtoCA.append(distance)
                    NtoCAtotal += distance
                elif curr.getName() == 'CA':
                    CAtoC.append(distance)
                    CAtoCtotal += distance
                    #get distance to adjacent CA atom
                    CA2 = curr.getCoords()  #current CA atom coordinates
                    if len(CA1) > 0:
                        CAtoCAdistance = math.sqrt((CA2[0] - CA1[0]) ** 2 + (CA2[1] - CA1[1]) ** 2 + (CA2[2] - CA1[2]) ** 2)    #distance between current and previous CA atoms
                        CAtoCA.append(CAtoCAdistance)
                        CAtoCAtotal += CAtoCAdistance
                    CA1 = CA2.copy()

                elif curr.getName() == 'C':
                    CtoN.append(distance)
                    CtoNtotal += distance
                if curr is not self.head:
                    dotP = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]    #dot product to find bond angle
                    magv1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)   #get magnitudes of both vectors
                    magv2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                    angle = math.acos(dotP/(magv1 * magv2))     #get angles (in radians)

                    norm2 = [v1[1] * v2[2] - v1[2] * v2[1], -1 * (v1[0] * v2[2] - v1[2] * v2[0]), v1[0] * v2[1] - v1[1] * v2[0]]  # current cross product of 2 vectors which will be used to find torsional angles
                    if len(norm1) > 0 and (curr.getRes() == 29 or curr.getRes() == 30):
                        dotPNorms = norm1[0] * norm2[0] + norm1[1] * norm2[1] + norm1[2] * norm2[2]     #dot product of vectors normal to the planes
                        magNorm1 = math.sqrt(norm1[0]**2 + norm1[1]**2 + norm1[2]**2)   #magnitude of first normal vector
                        magNorm2 = math.sqrt(norm2[0]**2 + norm2[1]**2 + norm2[2]**2)   #magnitude of second normal vector
                        torAngle = math.acos(dotPNorms/(magNorm1 * magNorm2))   #torsional angle

                    if curr.getName() == 'N':
                        angsCtoNtoCA.append(angle)
                        angtotalCtoNtoCA += angle
                        if curr.getRes() == 30:
                            self.stats['omegaAngle'] = torAngle     #omega torional angle
                    elif curr.getName() == 'CA':
                        angsNtoCAtoC.append(angle)
                        angtotalNtoCAtoC += angle
                        if curr.getRes() == 30:
                            self.stats['phiAngle'] = torAngle
                    elif curr.getName() == 'C':
                        angsCAtoCtoN.append(angle)
                        angtotalCAtoCtoN += angle
                        if curr.getRes() == 30:
                            self.stats['psiAngle'] = torAngle

                    norm1 = norm2.copy()

                v1 = v2.copy()
            curr = curr.nxt

        #get mean
        self.stats['meanNtoCA'] = NtoCAtotal/len(NtoCA)     #mean bond length of N to CA
        self.stats['meanCAtoC'] = CAtoCtotal/len(CAtoC)     #mean bond length of CA to C
        self.stats['meanCtoN'] = CtoNtotal/len(CtoN)    #mean bond length of C to N
        self.stats['meanAngleNCAC'] = angtotalNtoCAtoC/len(angsNtoCAtoC)    #mean bond angle of N to CA to C
        self.stats['meanAngleCACN'] = angtotalCAtoCtoN/len(angsCAtoCtoN)    #mean bond angle of CA to C to N
        self.stats['meanAngleCNCA'] = angtotalCtoNtoCA/len(angsCtoNtoCA)    #mean bond angle of C to N to CA
        self.stats['meanCAtoCA'] = CAtoCAtotal/len(CAtoCA)  #mean distance between adjacent CA atoms

        #get standard deviations
        varianceNtoCA = 0
        for lengths in NtoCA:
            varianceNtoCA += (lengths - self.stats['meanNtoCA']) ** 2
        self.stats['stdNtoCA'] = math.sqrt(varianceNtoCA/len(NtoCA))

        varianceCAtoC = 0
        for lengths in CAtoC:
            varianceCAtoC += (lengths - self.stats['meanCAtoC']) ** 2
        self.stats['stdCAtoC'] = math.sqrt(varianceCAtoC/len(CAtoC))

        varianceCtoN = 0
        for lengths in CtoN:
            varianceCtoN += (lengths - self.stats['meanCtoN']) ** 2
        self.stats['stdCtoN'] = math.sqrt(varianceCtoN/len(CtoN))

        varianceangNCAC = 0
        for ang in angsNtoCAtoC:
            varianceangNCAC += (ang - self.stats['meanAngleNCAC']) ** 2
        self.stats['stdAngleNCAC'] = math.sqrt(varianceangNCAC/len(angsNtoCAtoC))

        varianceangCACN = 0
        for ang in angsCAtoCtoN:
            varianceangCACN += (ang - self.stats['meanAngleCACN']) ** 2
        self.stats['stdAngleCACN'] = math.sqrt(varianceangCACN/len(angsCAtoCtoN))

        varianceangCNCA = 0
        for ang in angsCtoNtoCA:
            varianceangCNCA += (ang - self.stats['meanAngleCNCA']) ** 2
        self.stats['stdAngleCNCA'] = math.sqrt(varianceangCNCA/len(angsCtoNtoCA))

        varianceCAtoCA = 0
        for CAtoCAdistance in CAtoCA:
            varianceCAtoCA += (CAtoCAdistance - self.stats['meanCAtoCA']) ** 2
        self.stats['stdCAtoCA'] = math.sqrt(varianceCAtoCA/len(CAtoCA))


    def printList(self):
        curr = self.head
        while curr.getName() is not None:
            print('Name: ' + curr.getName())
            print('Coordinates (x,y,z): ' + str(curr.getCoords()))
            print()
            curr = curr.nxt
        print('Mean and standard deviation bond lengths:')
        print('Mean N-CA: ' + str(self.stats['meanNtoCA']))
        print('std N-CA: ' + str(self.stats['stdNtoCA']))
        print('Mean CA-C: ' + str(self.stats['meanCAtoC']))
        print('std CA-C: ' + str(self.stats['stdCAtoC']))
        print('Mean C-N: ' + str(self.stats['meanCtoN']))
        print('std C-N: ' + str(self.stats['stdCtoN']))
        print()
        print('Mean and std bond angles (radians):')
        print('Mean N-CA-C angle: ' + str(self.stats['meanAngleNCAC'])) #1.921
        print('std N-CA-C angle: ' + str(self.stats['stdAngleNCAC']))   #
        print('Mean CA-C-N angle: ' + str(self.stats['meanAngleCACN']))
        print('std CA-C-N angle: ' + str(self.stats['stdAngleCACN']))
        print('Mean C-N-CA angle: ' + str(self.stats['meanAngleCNCA']))
        print('std C-N-CA angle: ' + str(self.stats['stdAngleCNCA']))
        print()
        print('Mean distance between adjacent CA atoms: ' + str(self.stats['meanCAtoCA']))
        print('std distance between adjacent CA atoms: ' + str(self.stats['stdCAtoCA']))
        print()
        print('Torsional angles in residue 30 (in radians):')
        print('Phi angle: ' + str(self.stats['phiAngle']))
        print('Psi angle: ' + str(self.stats['psiAngle']))
        print('Omega angle: ' + str(self.stats['omegaAngle']))


if __name__ == '__main__':
    atoms = AtomList()
    atoms.readPDB('2GB1')
    atoms.getStats()
    atoms.printList()
