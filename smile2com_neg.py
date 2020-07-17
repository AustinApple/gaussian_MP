from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import openbabel
import pandas as pd
import os
import argparse
import numpy as np

list_smi = ['Cn1on1PC#N', 'FC(F)(F)N1CO1', 'CC1(C(F)(F)F)CO1', 'CC1(C(F)(F)F)OC1F', 'CC1(C(F)(F)F)OC1I', 
            'CC1(COF)CO1', 'CC1(F)CO1', 'CC1(F)OC1F', 'CC1(F)OCO1', 'CC1(IF)CO1', 'CS(=O)(=O)C(C1CO1)C(F)(F)F', 
            'FC(F)(F)C(C1CO1)C(F)(F)F', 'FC(F)(F)C1(CI)CO1', 'FC(F)(F)C1(F)CCO1', 'FC(F)(F)C1OCCO1', 
            'FC(F)(F)CC1(C(F)(F)F)CO1', 'FC(F)(F)CC1(F)CO1', 'FC(F)(F)CC1(I)CO1', 'FC(F)C1CO1', 'FC1(F)CCO1', 
            'FC1OCC1C(F)(F)F', 'FC1OCCO1', 'COC(=O)OC1CB1C(F)(F)F', 'CC(C)C1COC(=O)O1', 'FC(F)(F)OC1OI1', 
            'CC1(C(F)(F)F)CC(C(F)(F)F)C1', 'CC1(C(F)(F)F)CC1CC(F)(F)F', 'CC1(F)COC(=O)O1', 'CC1(I)COC(=O)O1', 
            'CC1CC(C(F)(F)F)C1F', 'CC1COC(=O)O1', 'CC1OC(=O)O1', 'FC(F)(F)C1CIO1', 'FC(F)(F)C1CO1', 
            'FC(F)(F)CC1CC(C(F)(F)F)C1', 'FC(F)(F)CCC1(C(F)(F)F)CC1', 'FC(F)(F)IC1CO1', 'FC1(F)CO1', 
            'FC1CC1CCC(F)(F)F', 'FCC1CC(C(F)(F)F)C1', 'FCC1CC1CC(F)(F)F', 'O=C(NBC1CO1)C(F)(F)F', 
            'CC(C)(C)n1on1F', 'CC(C)C1(F)COC(=O)O1', 'CC(C)Cn1on1F', 'CC(CF)C1COC(=O)O1', 'CC(CI)C1COC(=O)O1', 
            'CC(CI)C1OC(=O)O1', 'CC1(C(F)(F)F)CC1OC#N', 'CC1(C(F)F)CO1', 'CC1(C)OC(=O)OC1C(F)(F)F', 
            'CC1(C)OC1(F)C(F)(F)F', 'CC1(CC(F)(F)F)COC(=O)O1', 'CC1(CF)OC1F', 'CC1(CI)COC(=O)O1', 
            'CC1(COC(F)(F)F)COC(=O)C1', 'CC1C(F)OC1C(F)(F)F', 'CC1C(F)OC1F', 'CC1CC1(F)OC(F)(F)F', 
            'FC(F)(F)CC1(CC(F)(F)F)CO1', 'FC(F)(F)CC1(OC(F)(F)F)CC1', 'FC(F)(F)OC(C1CC1)C(F)(F)F', 
            'FC(F)(F)OC1(I)CO1', 'FC(F)(F)OCC1(C(F)(F)F)CC1', 'FC(F)(F)OCC1CC1C(F)(F)F', 'FC1(F)OCO1', 
            'FC1OCC1CC(F)(F)F', 'FC1OCO1', 'FCC1(CF)CO1', 'FIC1OCO1', 'FOC1OCCCO1', 'CC(C)(F)CC1OC(=O)O1', 
            'CC1(C(F)(F)F)CC(F)C1', 'CC1(C(F)(F)F)CCC1C(F)(F)F', 'CC1C(C(F)(F)F)CC1C(F)(F)F', 
            'CCC1(C(F)(F)F)CC1C(F)(F)F', 'CCC1CC1(F)C(F)(F)F', 'FC(F)(F)CC1(C(F)(F)F)CCC1', 'FC1CCC1CC(F)(F)F', 
            'FCC1(C(F)(F)F)CCC1', 'o1oo1', 'CC(C)COC(=O)n1oo1', 'CC(OF)n1oo1', 'CC1(C(F)(F)F)CIO1', 
            'FC(F)(F)CIC1CO1', 'FC(F)(F)COn1oo1', 'FC(F)(F)N1CIO1', 'CC(F)C1COC(=O)O1', 'CC(I)C1COC(=O)O1', 
            'CC1(C)OC(=O)O1', 'CC1(C)OC(=O)OI1', 'CC1(CI)OC(=O)O1', 'CCC1(F)COC(=O)O1', 'CCC1(F)OC(=O)O1', 
            'CCC1(I)COC(=O)O1', 'CCC1COC(=O)O1', 'CCC1OC(=O)O1', 'O=C1OC(CCI)O1', 'Fp1oo1']


for task_id, smi in enumerate(list_smi):
    print(smi)
    m = Chem.MolFromSmiles(smi)     
    m2=Chem.AddHs(m)
    AllChem.EmbedMolecule(m2)
    AllChem.MMFFOptimizeMolecule(m2)

    print(Chem.MolToMolBlock(m2),file=open(str(task_id)+"_neg.mol",'w+'))

    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("mol", "com")
    mol = openbabel.OBMol()

    obconversion.ReadFile(mol,str(task_id)+"_neg.mol")
    obconversion.WriteFile(mol,str(task_id)+"_neg.com")

    os.remove(str(task_id)+"_neg.mol")

    with open(str(task_id)+"_neg.com" ,'r+') as f_com:
        text = f_com.readlines()[5:]
    with open(str(task_id)+"_neg.com" ,'w+') as f_com:
        f_com.write("%mem=256MW\n# Opt Freq B3LYP/6-31+G* \
                     \n\n"+smi+"\n\n-1 2\n")
    with open(str(task_id)+"_neg.com" ,'a+') as f_com:
        for line in text:
            f_com.write(line)
        # f_com.write("eps=78.3553\n\n")
