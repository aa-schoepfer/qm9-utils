"""
QM9
    .na                 #number of atoms (int)
    .tag                #string (str)
    .i                  #index, integer identifier (int)
    .A                  #Rotaional Constant (float)
    .B                  #Rotaional Constant (float)
    .C                  #Rotaional Constant (float)
    .mu                 #Dipole moment (float)
    .alpha              #Isotropic polarizability (float)
    .homo               #Energy of HOMO (float)
    .lumo               #Energy of LUMO (float)
    .gap                #Gap (eLUMO - eHOMO) (float)
    .ese                #Electronic spacial extent (float)
    .zpve               #Zero point vibrational energy (float)
    .U0                 #Internal energy at 0K (float)
    .U                  #Internal energy at 298.15K (float)
    .H                  #Enthalpy at 298.15K (float)
    .G                  #Free energy at 298.15K (float)
    .Cv                 #Heat capacity at 298.15K (float)
    .mols               #Mol object (obj)
    [i]
        .atom           #Element type (str)
        .x              #X coordinate (float)
        .y              #Y coordinate (float)
        .z              #Z coordinate (float)
        .mulliken       #Muliken partial charge (float)
    .ncs                #Number of chrages (np.array)
    .coordinates        #Coordinates (np.array)
    .mullikens          #Muliken partial charges (np.array)
    .atomtypes          #Dictionary with frequency per atomtype (dict)
    .mw                 #Aproximative molar weight (float)
    .hvfs               #Harmonic vibrational frequencies (np.array)
    [i]
    .smiles             #SMILES
        .gdb17          #(str)
        .b3lyp          #(str)
    .inchis             #InChi
        .corina         #(str)
        .b3lyp          #(str)
    .representation     #To use with QML code

"""

import numpy as np

nuclear_charge = {
 'H'  :     1,
 'He' :     2,
 'Li' :     3,
 'Be' :     4,
 'B'  :     5,
 'C'  :     6,
 'N'  :     7,
 'O'  :     8,
 'F'  :     9,
 'Ne' :    10,
 'Na' :    11,
 'Mg' :    12,
 'Al' :    13,
 'Si' :    14,
 'P'  :    15,
 'S'  :    16,
 'Cl' :    17,
 'Ar' :    18,
 'K'  :    19,
 'Ca' :    20,
 'Sc' :    21,
 'Ti' :    22,
 'V'  :    23,
 'Cr' :    24,
 'Mn' :    25,
 'Fe' :    26,
 'Co' :    27,
 'Ni' :    28,
 'Cu' :    29,
 'Zn' :    30,
 'Ga' :    31,
 'Ge' :    32,
 'As' :    33,
 'Se' :    34,
 'Br' :    35,
 'Kr' :    36,
 'Rb' :    37,
 'Sr' :    38,
 'Y'  :    39,
 'Zr' :    40,
 'Nb' :    41,
 'Mo' :    42,
 'Tc' :    43,
 'Ru' :    44,
 'Rh' :    45,
 'Pd' :    46,
 'Ag' :    47,
 'Cd' :    48,
 'In' :    49,
 'Sn' :    50,
 'Sb' :    51,
 'Te' :    52,
 'I'  :    53,
 'Xe' :    54,
 'Cs' :    55,
 'Ba' :    56,
 'La' :    57,
 'Ce' :    58,
 'Pr' :    59,
 'Nd' :    60,
 'Pm' :    61,
 'Sm' :    62,
 'Eu' :    63,
 'Gd' :    64,
 'Tb' :    65,
 'Dy' :    66,
 'Ho' :    67,
 'Er' :    68,
 'Tm' :    69,
 'Yb' :    70,
 'Lu' :    71,
 'Hf' :    72,
 'Ta' :    73,
 'W'  :    74,
 'Re' :    75,
 'Os' :    76,
 'Ir' :    77,
 'Pt' :    78,
 'Au' :    79,
 'Hg' :    80,
 'Tl' :    81,
 'Pb' :    82,
 'Bi' :    83,
 'Po' :    84,
 'At' :    85,
 'Rn' :    86,
 'Fr' :    87,
 'Ra' :    88,
 'Ac' :    89,
 'Th' :    90,
 'Pa' :    91,
 'U'  :    92,
 'Np' :    93,
 'Pu' :    94,
 'Am' :    95,
 'Cm' :    96,
 'Bk' :    97,
 'Cf' :    98,
 'Es' :    99,
 'Fm' :   100,
 'Md' :   101,
 'No' :   102,
 'Lr' :   103,
 'Rf' :   104,
 'Db' :   105,
 'Sg' :   106,
 'Bh' :   107,
 'Hs' :   108,
 'Mt' :   109,
 'Ds' :   110,
 'Rg' :   111,
 'Cn' :   112,
 'Uuq':   114,
 'Uuh':   116}

class Mol:
    def __init__(self, mol):
        self.atom = mol[0]
        self.x = float(mol[1].replace("*^","e"))
        self.y = float(mol[2].replace("*^","e"))
        self.z = float(mol[3].replace("*^","e"))
        self.mulliken = float(mol[4].replace("*^","e"))

class Smile:
    def __init__(self, smiles):
        self.gdb17 = smiles[0]
        self.b3lyp = smiles[1]

class Inchi:
    def __init__(self, inchis):
        self.corina = inchis[0]
        self.b3lyp = inchis[1]

class QM9:

    def __init__(self, file):

        f = open(file)
        data = f.readlines()
        f.close()
        len_data = len(data)

        for index, line in enumerate(data):
            data[index] = line.split()
        self.na = int(data[0][0])
        self.tag = data[1][0]
        self.i = int(data[1][1])
        self.A = float(data[1][2])
        self.B = float(data[1][3])
        self.C = float(data[1][4])
        self.mu = float(data[1][5])
        self.alpha = float(data[1][6])
        self.homo = float(data[1][7])
        self.lumo = float(data[1][8])
        self.gap = float(data[1][9])
        self.ese = float(data[1][10])
        self.zpve = float(data[1][11])
        self.U0 = float(data[1][12])
        self.U = float(data[1][13])
        self.H = float(data[1][14])
        self.G = float(data[1][15])
        self.Cv = float(data[1][16])
        self.mols = data[2:len_data-3]
        self.ncs = np.empty(self.na, dtype=int)
        self.coordinates =  np.empty((self.na, 3), dtype=float)
        self.mullikens = np.empty(self.na, dtype=float)
        for index, mol in enumerate(self.mols):
            for i in range(1,5):
                mol[i] = mol[i].replace("*^","e")
            self.mols[index] = Mol(mol)
            self.ncs[index] = nuclear_charge[mol[0]]
            self.coordinates[index] = np.asarray(mol[1:4], dtype=float)
            self.mullikens[index] = float(mol[4])
        self.atomtypes = {}
        for mol in self.mols:
            if mol.atom in self.atomtypes:
                self.atomtypes[mol.atom] += 1
            else:
                self.atomtypes[mol.atom] = 1
        self.mw = 0
        am = {'H': 1.007947, 'C': 12.0111, 'O':15.99943, 'N':14.006747, 'F':18.99840329}
        for mol in self.mols:
            self.mw += am[mol.atom]
        self.hvfs = data[len_data-3]
        self.hvfs = np.array([float(hvf) for hvf in self.hvfs])
        self.smiles = Smile(data[len_data-2])
        self.inchis = Inchi(data[len_data-1])
        self.representation = np.asarray([], dtype = float)
