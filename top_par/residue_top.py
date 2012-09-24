### Alias for residues

residue={
## residues
'ALA'     : 'reALA'    ,
'ARG'     : 'reARG'    ,
'ASN'     : 'reASN'    ,
'ASP'     : 'reASP'    ,
'CYS'     : 'reCYS'    ,
'GLU'     : 'reGLU'    ,
'GLN'     : 'reGLN'    ,
'GLY'     : 'reGLY'    ,
'HIS'     : 'reHIS'    ,
'HSD'     : 'reHIS'    ,
'ILE'     : 'reILE'    ,
'LEU'     : 'reLEU'    ,
'LYS'     : 'reLYS'    ,
'LYSH'    : 'reLYS'    ,
'MET'     : 'reMET'    ,
'PHE'     : 'rePHE'    ,
'PRO'     : 'rePRO'    ,
'SER'     : 'reSER'    ,
'THR'     : 'reTHR'    ,
'TRP'     : 'reTRP'    ,
'TYR'     : 'reTYR'    ,
'VAL'     : 'reVAL'    ,
# terminals
'ACE'     : 'reACE'    ,
'NME'     : 'reNME'    ,
'NAC'     : 'reNME'    ,
'NHE'     : 'reNHE'    ,
'NH2'     : 'reNHE'    ,
# water
'SOL'     : 'reSOL3'    ,
'HOH'     : 'reSOL3'    ,
'TIP'     : 'reSOL3'    ,
'TIP3'    : 'reSOL3'    ,
'HO4'     : 'reSOL4'    ,
'TIP4'    : 'reSOL4'    ,
'HO5'     : 'reSOL5'    ,
'TIP5'    : 'reSOL5'    ,
'SWM'     : 'reSOL5'    ,
# ions
'NA'      : 'reNA'      ,
'K'       : 'reK'       ,
'LI'      : 'reLI'      ,
'CL'      : 'reCL'      
}

#residue_atoms={'amino':[A,B,C],...}
#covalent_bonds={'amino':[[A,B],[A,C]],...}

residue_atoms={}
covalent_bonds={}
terminal_atoms={}
terminal_bonds={}

###############################
###############################

############ PEPTIDES::

######## Aminoacids

### ALA:

residue_atoms['reALA']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atHB3',
'atC',
'atO'
]
 
covalent_bonds['reALA']=[
['atN'   ,'atH'    ],
['atN'   ,'atCA'   ],
['atCA'  ,'atHA'   ],
['atCA'  ,'atCB'   ],
['atCA'  ,'atC'    ],
['atCB'  ,'atHB1'  ],
['atCB'  ,'atHB2'  ],
['atCB'  ,'atHB3'  ],
['atC'   ,'atO'    ] 
]

### ARG:

residue_atoms['reARG']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atHD1', 
'atHD2', 
'atNE',
'atHE',
'atCZ',
'atNH1',  
'atHH11', 
'atHH12', 
'atNH2',  
'atHH21', 
'atHH22', 
'atC',
'atO'
]

covalent_bonds['reARG']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atHG1'  ], 
['atCG'  ,'atHG2'  ], 
['atCG'  ,'atCD'   ], 
['atCD'  ,'atHD1'  ], 
['atCD'  ,'atHD2'  ], 
['atCD'  ,'atNE'   ], 
['atNE'  ,'atHE'   ], 
['atNE'  ,'atCZ'   ], 
['atCZ'  ,'atNH1'  ], 
['atCZ'  ,'atNH2'  ], 
['atNH1' ,'atHH11' ], 
['atNH1' ,'atHH12' ], 
['atNH2' ,'atHH21' ], 
['atNH2' ,'atHH22' ], 
['atC'   ,'atO'    ] 
]

### ASN:

residue_atoms['reASN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1',  
'atND2',  
'atHD21', 
'atHD22', 
'atC',
'atO'
]

covalent_bonds['reASN']=[
['atN'   ,'atH'     ], 
['atN'   ,'atCA'    ], 
['atCA'  ,'atHA'    ], 
['atCA'  ,'atCB'    ], 
['atCA'  ,'atC'     ], 
['atCB'  ,'atHB1'   ], 
['atCB'  ,'atHB2'   ], 
['atCB'  ,'atCG'    ], 
['atCG'  ,'atOD1'   ], 
['atCG'  ,'atND2'   ], 
['atND2' ,'atHD21'  ], 
['atND2' ,'atHD22'  ], 
['atC'   ,'atO'     ] 
]

### ASP:

residue_atoms['reASP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atOD1', 
'atOD2', 
'atC',
'atO'
]

covalent_bonds['reASP']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atOD1'  ], 
['atCG'  ,'atOD2'  ], 
['atC'   ,'atO'    ] 
]

### CYS and CYSH [HG1]:

residue_atoms['reCYS']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atSG',
'atHG1',
'atC',
'atO'
]
 
covalent_bonds['reCYS']=[
['atN'   ,'atH'    ],
['atN'   ,'atCA'   ],
['atCA'  ,'atHA'   ],
['atCA'  ,'atCB'   ],
['atCA'  ,'atC'    ],
['atCB'  ,'atHB1'  ],
['atCB'  ,'atHB2'  ],
['atCB'  ,'atSG'   ],
['atSG'  ,'atHG1'  ],
['atC'   ,'atO'    ] 
]

### GLU:

residue_atoms['reGLU']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atOE1', 
'atOE2', 
'atC',
'atO'
]

covalent_bonds['reGLU']=[
['atN'   , 'atH'   ], 
['atN'   , 'atCA'  ], 
['atCA'  , 'atHA'  ], 
['atCA'  , 'atCB'  ], 
['atCA'  , 'atC'   ], 
['atCB'  , 'atHB1' ], 
['atCB'  , 'atHB2' ], 
['atCB'  , 'atCG'  ], 
['atCG'  , 'atHG1' ], 
['atCG'  , 'atHG2' ], 
['atCG'  , 'atCD'  ], 
['atCD'  , 'atOE1' ], 
['atCD'  , 'atOE2' ], 
['atC'   , 'atO'   ] 
]

### GLN:

residue_atoms['reGLN']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atHG1', 
'atHG2', 
'atCD',
'atOE1', 
'atNE2', 
'atHE21',
'atHE22',
'atC',
'atO'
]

covalent_bonds['reGLN']=[
['atN'   , 'atH'   ], 
['atN'   , 'atCA'  ], 
['atCA'  , 'atHA'  ], 
['atCA'  , 'atCB'  ], 
['atCA'  , 'atC'   ], 
['atCB'  , 'atHB1' ], 
['atCB'  , 'atHB2' ], 
['atCB'  , 'atCG'  ], 
['atCG'  , 'atHG1' ], 
['atCG'  , 'atHG2' ], 
['atCG'  , 'atCD'  ], 
['atCD'  , 'atOE1' ], 
['atCD'  , 'atNE2' ], 
['atNE2' , 'atHE21'], 
['atNE2' , 'atHE22'], 
['atC'   , 'atO'   ] 
]

### GLY:

residue_atoms['reGLY']=[
'atN',
'atH',
'atCA',
'atHA1',
'atHA2',
'atC',
'atO'
]

covalent_bonds['reGLY']=[
['atN'   , 'atH'   ],
['atN'   , 'atCA'  ],
['atCA'  , 'atHA1' ],
['atCA'  , 'atHA2' ],
['atCA'  , 'atC'   ],
['atC'   , 'atO'   ] 
]

### HISE (ND1 no H, NE2 with H),
### HISD (ND1 with H, NE2 no H),
### HISH (ND1 with H, NE2 with H),
### All included in HIS:

residue_atoms['reHIS']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2',
'atCG',
'atND1', 
'atHD1',
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atNE2', 
'atHE2', 
'atC',
'atO'
]

covalent_bonds['reHIS']=[
['atN'   ,'atH'    ], 
['atN'   ,'atCA'   ], 
['atCA'  ,'atHA'   ], 
['atCA'  ,'atCB'   ], 
['atCA'  ,'atC'    ], 
['atCB'  ,'atHB1'  ], 
['atCB'  ,'atHB2'  ], 
['atCB'  ,'atCG'   ], 
['atCG'  ,'atND1'  ], 
['atCG'  ,'atCD2'  ], 
['atND1' ,'atHD1'  ], 
['atND1' ,'atCE1'  ], 
['atCD2' ,'atHD2'  ], 
['atCD2' ,'atNE2'  ], 
['atCE1' ,'atHE1'  ], 
['atCE1' ,'atNE2'  ], 
['atNE2' ,'atHE2'  ], 
['atC'   ,'atO'    ] 
]

### ILE:

residue_atoms['reILE']=[
'atN', 
'atH', 
'atCA', 
'atHA', 
'atCB', 
'atHB', 
'atCG1', 
'atHG11', 
'atHG12', 
'atCG2', 
'atHG21', 
'atHG22', 
'atHG23', 
'atCD', 
'atHD1', 
'atHD2', 
'atHD3', 
'atC', 
'atO' 
]

covalent_bonds['reILE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB'   ], 
['atCB',  'atCG1'  ], 
['atCB',  'atCG2'  ], 
['atCG1', 'atHG11' ], 
['atCG1', 'atHG12' ], 
['atCG1', 'atCD'   ], 
['atCG2', 'atHG21' ], 
['atCG2', 'atHG22' ], 
['atCG2', 'atHG23' ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atHD3'  ], 
['atC',   'atO'    ] 
]

### LEU:

residue_atoms['reLEU']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1',
'atHB2',
'atCG',
'atHG1',
'atCD1',
'atHD11', 
'atHD12', 
'atHD13', 
'atCD2',
'atHD21', 
'atHD22', 
'atHD23', 
'atC',
'atO'
]

covalent_bonds['reLEU']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD11' ], 
['atCD1', 'atHD12' ], 
['atCD1', 'atHD13' ], 
['atCD2', 'atHD21' ], 
['atCD2', 'atHD22' ], 
['atCD2', 'atHD23' ], 
['atC',   'atO'    ] 
]

### LYS with LYSH (HZ3):

residue_atoms['reLYS']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atCD', 
'atHD1',
'atHD2',
'atCE', 
'atHE1',
'atHE2',
'atNZ', 
'atHZ1',
'atHZ2',
'atHZ3',
'atC',
'atO'
]

covalent_bonds['reLYS']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atHG1' ], 
['atCG',  'atHG2' ], 
['atCG',  'atCD'  ], 
['atCD',  'atHD1' ], 
['atCD',  'atHD2' ], 
['atCD',  'atCE'  ], 
['atCE',  'atHE1' ], 
['atCE',  'atHE2' ], 
['atCE',  'atNZ'  ], 
['atNZ',  'atHZ1' ], 
['atNZ',  'atHZ2' ], 
['atNZ',  'atHZ3' ], 
['atC',   'atO'   ] 
]

### MET:

residue_atoms['reMET']=[
'atN',
'atH', 
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atSD', 
'atCE', 
'atHE1',
'atHE2',
'atHE3',
'atC',
'atO'
]

covalent_bonds['reMET']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atSD'   ], 
['atSD',  'atCE'   ], 
['atCE',  'atHE1'  ], 
['atCE',  'atHE2'  ], 
['atCE',  'atHE3'  ], 
['atC',   'atO'    ] 
]

### PHE:

residue_atoms['rePHE']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD1', 
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ',
'atHZ',
'atC',
'atO'
]

covalent_bonds['rePHE']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD1'  ], 
['atCD1', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atHZ'   ], 
['atC',   'atO'    ] 
]

### PRO:

residue_atoms['rePRO']=[
'atN',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atHG1',
'atHG2',
'atCD', 
'atHD1',
'atHD2',
'atC', 
'atO'
]

covalent_bonds['rePRO']=[
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atHG1'  ], 
['atCG',  'atHG2'  ], 
['atCG',  'atCD'   ], 
['atCD',  'atHD1'  ], 
['atCD',  'atHD2'  ], 
['atCD',  'atN'    ], 
['atC',   'atO'    ] 
]

### SER:

residue_atoms['reSER']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atOG', 
'atHG1', 
'atC',
'atO'
]

covalent_bonds['reSER']=[
['atN',  'atH'    ], 
['atN',  'atCA'   ], 
['atCA', 'atHA'   ], 
['atCA', 'atCB'   ], 
['atCA', 'atC'    ], 
['atCB', 'atHB1'  ], 
['atCB', 'atHB2'  ], 
['atCB', 'atOG'   ], 
['atOG', 'atHG1'  ], 
['atC',  'atO'    ] 
]

### THR:

residue_atoms['reTHR']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB',
'atOG1',
'atHG1',
'atCG2',
'atHG21',
'atHG22', 
'atHG23', 
'atC',
'atO'
]

covalent_bonds['reTHR']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atOG1'   ], 
['atCB',  'atCG2'   ], 
['atOG1', 'atHG1'   ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

### TRP:

residue_atoms['reTRP']=[
'atN',
'atH',
'atCA',
'atHA',
'atCB',
'atHB1', 
'atHB2', 
'atCG',
'atCD1', 
'atHD1', 
'atCD2', 
'atNE1', 
'atHE1', 
'atCE2', 
'atCE3', 
'atHE3', 
'atCZ2', 
'atHZ2', 
'atCZ3', 
'atHZ3', 
'atCH2', 
'atHH2', 
'atC',
'atO'
]

covalent_bonds['reTRP']=[
['atN',   'atH'   ], 
['atN',   'atCA'  ], 
['atCA',  'atHA'  ], 
['atCA',  'atCB'  ], 
['atCA',  'atC'   ], 
['atCB',  'atHB1' ], 
['atCB',  'atHB2' ], 
['atCB',  'atCG'  ], 
['atCG',  'atCD1' ], 
['atCG',  'atCD2' ], 
['atCD1', 'atHD1' ], 
['atCD1', 'atNE1' ], 
['atCD2', 'atCE2' ], 
['atCD2', 'atCE3' ], 
['atNE1', 'atHE1' ], 
['atNE1', 'atCE2' ], 
['atCE2', 'atCZ2' ], 
['atCE3', 'atHE3' ], 
['atCE3', 'atCZ3' ], 
['atCZ2', 'atHZ2' ], 
['atCZ2', 'atCH2' ], 
['atCZ3', 'atHZ3' ], 
['atCZ3', 'atCH2' ], 
['atCH2', 'atHH2' ], 
['atC',   'atO'   ] 
]

### TYR:

residue_atoms['reTYR']=[
'atN',
'atH',
'atCA', 
'atHA', 
'atCB', 
'atHB1',
'atHB2',
'atCG', 
'atCD1',
'atHD1', 
'atCD2', 
'atHD2', 
'atCE1', 
'atHE1', 
'atCE2', 
'atHE2', 
'atCZ', 
'atOH', 
'atHH', 
'atC',
'atO'
]

covalent_bonds['reTYR']=[
['atN',   'atH'    ], 
['atN',   'atCA'   ], 
['atCA',  'atHA'   ], 
['atCA',  'atCB'   ], 
['atCA',  'atC'    ], 
['atCB',  'atHB1'  ], 
['atCB',  'atHB2'  ], 
['atCB',  'atCG'   ], 
['atCG',  'atCD1'  ], 
['atCG',  'atCD2'  ], 
['atCD1', 'atHD1'  ], 
['atCD1', 'atCE1'  ], 
['atCD2', 'atHD2'  ], 
['atCD2', 'atCE2'  ], 
['atCE1', 'atHE1'  ], 
['atCE1', 'atCZ'   ], 
['atCE2', 'atHE2'  ], 
['atCE2', 'atCZ'   ], 
['atCZ',  'atOH'   ], 
['atOH',  'atHH'   ], 
['atC',   'atO'    ] 
]

### VAL:

residue_atoms['reVAL']=[
'atN',
'atH', 
'atCA',
'atHA',
'atCB',
'atHB',
'atCG1', 
'atHG11',
'atHG12',
'atHG13',
'atCG2', 
'atHG21',
'atHG22',
'atHG23',
'atC',
'atO'
]

covalent_bonds['reVAL']=[
['atN',   'atH'     ], 
['atN',   'atCA'    ], 
['atCA',  'atHA'    ], 
['atCA',  'atCB'    ], 
['atCA',  'atC'     ], 
['atCB',  'atHB'    ], 
['atCB',  'atCG1'   ], 
['atCB',  'atCG2'   ], 
['atCG1', 'atHG11'  ], 
['atCG1', 'atHG12'  ], 
['atCG1', 'atHG13'  ], 
['atCG2', 'atHG21'  ], 
['atCG2', 'atHG22'  ], 
['atCG2', 'atHG23'  ], 
['atC',   'atO'     ] 
]

######## Terminals

### ACE:

residue_atoms['reACE']=[
'atCH3',
'atHH31',
'atHH32',
'atHH33',
'atC',
'atO'
]

covalent_bonds['reACE']=[
['atCH3' , 'atHH31' ],
['atCH3' , 'atHH32' ],
['atCH3' , 'atHH33' ],
['atCH3' , 'atC'    ],
['atC'   , 'atO'    ] 
]

### NME:

residue_atoms['reNME']=[
'atN',
'atH',
'atCH3',
'atHH31',
'atHH32',
'atHH33'
]

covalent_bonds['reNME']=[
['atN'   ,'atH'    ],
['atN'   ,'atCH3'  ],
['atCH3' ,'atHH31' ],
['atCH3' ,'atHH32' ],
['atCH3' ,'atHH33' ] 
]

### NHE:

residue_atoms['reNHE']=[
'atN',
'atH1',
'atH2'
]

covalent_bonds['reNHE']=[
['atN'   ,'atH1'   ],
['atN'   ,'atH2'   ]
]


##### WATER:

### SOL3:

residue_atoms['reSOL3']=[
'atOW',
'atHW1',
'atHW2'
]

covalent_bonds['reSOL3']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### SOL4:

residue_atoms['reSOL4']=[
'atOW',
'atHW1',
'atHW2',
'atvir'
]

covalent_bonds['reSOL4']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]

### SOL5:

residue_atoms['reSOL5']=[
'atOW',
'atHW1',
'atHW2',
'atvir',
'atvir'
]

covalent_bonds['reSOL5']=[
['atOW'  ,'atHW1' ],
['atOW'  ,'atHW2' ]
]


##### IONS:

## NA:
residue_atoms['reNA']=[
'itNA'
]

covalent_bonds['reNA']=[
]

## K:
residue_atoms['reK']=[
'itK'
]

covalent_bonds['reK']=[
]

## LI:
residue_atoms['reLI']=[
'itLI'
]

covalent_bonds['reLI']=[
]

## CL:
residue_atoms['reCL']=[
'itCL'
]

covalent_bonds['reCL']=[
]


