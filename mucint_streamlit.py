from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdFingerprintGenerator, MACCSkeys, Descriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataStructs import cDataStructs
from io import StringIO
from mordred import Calculator, descriptors
import numpy as np
import pandas as pd
import seaborn as sns
import sys, os, shutil
import matplotlib.pyplot as plt
import streamlit as st
from streamlit_ketcher import st_ketcher
import time
import subprocess
from PIL import Image
import uuid
from filelock import Timeout, FileLock

def cooling_highlight(val):
   color = "red" if val < 50 else "green"                    
   return f'background-color: {color}'

calc = Calculator(descriptors, ignore_3D=False)

def fingerprint_rdk5(self) -> np.ndarray:
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=5,fpSize=16384)
    return fp_gen.GetCountFingerprintAsNumPy(self).astype(int)

def fingerprint_rdk7(self) -> np.ndarray:
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7,fpSize=16384)
    return fp_gen.GetCountFingerprintAsNumPy(self).astype(int)

def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    clean_mol = rdMolStandardize.Cleanup(mol) 
    parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
    uncharger = rdMolStandardize.Uncharger()
    uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)    
    te = rdMolStandardize.TautomerEnumerator()
    taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)     
    taut_uncharged_parent_clean_mol_addH = Chem.AddHs(taut_uncharged_parent_clean_mol)
    Chem.SanitizeMol(taut_uncharged_parent_clean_mol_addH)
    return taut_uncharged_parent_clean_mol_addH

def cooling_highlight(val):
   color = 'green' if val == 8 else "green" if val == 7 else "green" if val == 6 else "green" if val == "X1" else "yellow" if val == 5 else "yellow" if val == 4  else "yellow" if val == 3 else "red" if val == 2 else "red" if val == 1 else "red" if val == 0 else "red" if val == "X0" else "grey" if val == "AD" else "white"                    
   return f'background-color: {color}'

st.title('MUC2 interaction probability model')

with st.form(key='my_form_to_submit'):
    with st.expander("More information"):
        
        st.caption(""":black[Background]""")
        st.caption("""mucint predicts interactions of drugs with mucin, based on classifications derived from ¹H-NMR measurements with MUC2. It is based on a logistic regression model, using five different MACCS keys as descriptors.""")
        
        st.caption("""The software is hosted at our [github page](https://github.com/juppifluppi/mucint), licensed under MIT.""")
 
        st.caption("""Version 0.1 (11.06.2024)""")
 
    SMI = st.text_input('Enter [SMILES code](https://pubchem.ncbi.nlm.nih.gov//edit3/index.html) of drug to load', '') 
    
#    on = st.toggle('Use drawn structure')
#    with st.expander("SMILES editor"):
#        drawer = st_ketcher()
#        st.caption("Click on Apply to save the drawn structure as input.")
#
#    if on:
#        SMI=drawer


    on3 = st.toggle('Perform batch calculation',key="16")    
    with st.expander("Batch settings"):
        col1, col2 = st.columns(2)
        with col1:
            st.write("Names of compounds")
            NAMESx = st.text_area(label="Input names of compounds separated by linebreaks",key="17")
        with col2:
            st.write("SMILES codes")
            SMILESx = st.text_area(label="Input SMILES of compounds separated by linebreaks",key="18")

    
    emoji = ''
    label = ' Predict'    
    submit_button = st.form_submit_button(label=f'{emoji} {label}')

if submit_button:
# try:
    lock = FileLock("lockfile.lock")
    with st.spinner('WAITING IN QUEUE ...'):
        try:
            lock.acquire(timeout=20)
        except Timeout:
            os.remove("lockfile.lock")
            lock = FileLock("lockfile.lock")
            lock.acquire()

    try:
        if on3 is False:
    
            for es in ["descriptors.csv","results.csv","results2.csv"]:
                try:
                    os.remove(es)
                except:
                    pass
    
            SMI=SMI
                                 
            mol = standardize(SMI)
            maccskeys = MACCSkeys.GenMACCSKeys(mol)            

            with open("descriptors.csv","a") as f:
                for o in range(0,len(maccskeys)):
                    f.write("maccs_"+str(o)+"\t")    
                f.write("\n")  
                
            with open("descriptors.csv","a") as f:
                for o in range(0,len(maccskeys)):
                    f.write(str(maccskeys[o])+"\t") 
                f.write("\n")
                                                                                     
            process3=subprocess.Popen(["Rscript", "predict.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            process3.communicate()
           
            df2 = pd.read_csv(r'results.csv')
            df3 = pd.read_csv(r'results2.csv')
           
            col1, col2 = st.columns(2)
                    
            with col1: 
                st.write("MUC2 interaction probabilty: "+str(int(df2.iloc[0, 0]*100))+" %")
                st.write("In presence of bile: "+str(int(df3.iloc[0, 0]*100))+" %")               
                            
            with col2:
                im = Draw.MolToImage(Chem.MolFromSmiles(SMI),fitImage=True)
                st.image(im)
            
            for es in ["descriptors.csv","results.csv","results2.csv"]:
                try:
                    os.remove(es)
                except:
                    pass
   
        if on3 is True:

           SMILESx=SMILESx.split('\n')           
           NAMESx=NAMESx.split('\n')

           for es in ["descriptors.csv"]:
                try:
                    os.remove(es)
                except:
                    pass

           for cx in range(0,len(SMILESx)):
               SMI=SMILESx[cx]
                          
               mol = standardize(SMI)
               maccskeys = MACCSkeys.GenMACCSKeys(mol)            
    
               if cx == 0:

                   with open("descriptors.csv","a") as f:
                       for o in range(0,len(maccskeys)):
                           f.write("maccs_"+str(o)+"\t")    
                       f.write("\n")  
                
               with open("descriptors.csv","a") as f:
                   for o in range(0,len(maccskeys)):
                       f.write(str(maccskeys[o])+"\t") 
                   f.write("\n")
            
           process3=subprocess.Popen(["Rscript", "predict.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
           process3.communicate()
                                           
           df2 = pd.read_csv(r'results.csv')
           df3 = pd.read_csv(r'results2.csv')        
        
           dfx = pd.DataFrame(columns=['Compound', "MUC2 interaction probability", "bile"])
           dfx["Compound"]=NAMESx
           dfx["MUC2 interaction probability"]=(df2.iloc[:, 0].astype(float))*100
           dfx["MUC2 interaction probability"]=dfx.iloc[:, 1].astype(int)

           dfx["bile"]=(df3.iloc[:, 0].astype(float))*100
           dfx["bile"]=dfx.iloc[:, 2].astype(int)
    
           #dfx.reset_index(inplace=True)               
           st.dataframe(dfx.style.applymap(cooling_highlight,subset=["MUC2 interaction probability", "bile"]))    
    
    finally:
        lock.release()

#    except:
#        st.write("Something went wrong. Cannot parse molecules! Please verify your structures.")  

