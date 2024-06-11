from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdFingerprintGenerator, Descriptors, Draw
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

with st.form(key='my_form_to_submit'):
    with st.expander("More information"):
        
        st.caption(""":black[Background]""")
        st.caption("""mucint predicts interactions of drugs with mucin, based on classifications""")
        
        st.caption("""Details can be found in our [publication](https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.4c00086) or the open-access [preprint](https://doi.org/10.26434/chemrxiv-2024-l5kvc) of it.""")
        
        st.caption("""The software is hosted at our [github page](https://github.com/juppifluppi/mucint), licensed under MIT.""")
 
        st.caption("""Version 0.1 (11.06.2024)""")
 
    SMI = st.text_input('Enter [SMILES code](https://pubchem.ncbi.nlm.nih.gov//edit3/index.html) of drug to load', '') 
    
    on = st.toggle('Use drawn structure',key="13")
    with st.expander("SMILES editor"):
        drawer = st_ketcher(key="12")
        st.caption("Click on Apply to save the drawn structure as input.")

    if on:
        SMI=drawer
    
    emoji = ''
    label = ' Predict'    
    submit_button = st.form_submit_button(label=f'{emoji} {label}')

if submit_button:
# try:
    lock = FileLock("lockfile.lock")
    with st.spinner('WAITING IN QUEUE ...'):
        try:
            lock.acquire(timeout=420)
        except Timeout:
            os.remove("lockfile.lock")
            lock = FileLock("lockfile.lock")
            lock.acquire()
    
    try:
        for es in ["descriptors.csv"]:
            try:
                os.remove(es)
            except:
                pass
    
        SMI=SMI
                                 
        mol = standardize(SMILES[molecule])
        maccskeys = MACCSkeys.GenMACCSKeys(mol)            

        with open("descriptors.csv","a") as f:
            for o in range(0,len(maccskeys)):
                f.write("maccs_"+str(o)+"\t")    
            f.write("\n")  
            
        with open("descriptors.csv","a") as f:
            for o in range(0,len(maccskeys)):
                f.write(str(maccskeys[o])+"\t") 
            f.write("\n")
                                                                                 
        process3 = subprocess.Popen(["Rscript", "predict.R"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        result3 = process3.communicate()
                                           
                    df2 = pd.read_csv(r'fin_results2.csv')
                    df2 = df2.rename(columns={0: "POL", 1: "DF", 2: "LC", 3: "LE"})
            
                    SDc = ((df2["DF"])*((df2["LE"])/100))
                    SDc2 = (((df2["LC"]/100)*(-1)*10)/((df2["LC"]/100)-1))
                    SDcx = ((SDc+SDc2)/2)
                    
                    if len(SMI2) > 2:
                        SDc = ((df2["DF"])*((df2["LE"])/100))
                        SDc2 = (((df2["LC"]/100)*(-1)*10)/((df2["LC"]/100)-1))
                        SDcx = ((SDc+SDc2)/2)
                    
            
                    calcLE=(SDcx/df2["DF"])*100
                    calcLC=(SDcx/(SDcx+10))*100
                    
                    df3={'POL' : df2["POL"], 'DF' : df2["DF"], 'SD': SDcx}
                    df3=pd.DataFrame(df3,columns=["POL","DF","SD"])
                    
                    df4={'POL' : df2["POL"], 'DF' : df2["DF"], 'SD': SDc}
                    df4=pd.DataFrame(df4,columns=["POL","DF","SD"])
            
                    df5={'POL' : df2["POL"], 'DF' : df2["DF"], 'SD': SDc2}
                    df5=pd.DataFrame(df5,columns=["POL","DF","SD"])
            
                    df6={'POL' : df2["POL"], 'DF' : df2["DF"], 'LE': calcLE,'LC': calcLC}
                    df6=pd.DataFrame(df6,columns=["POL","DF","LE","LC"])
                    
                    custom_palette = sns.color_palette("deep")
            
                    max_indexes = SDcx[SDcx == max(SDcx)].index.tolist()
                    
                    col1, col2 = st.columns(2)
            
                    finalLE=round((round(max(SDcx),1)/df3.loc[SDcx.idxmax(), "DF"])*100,0)
                    finalLC=round((round(max(SDcx),1)/(df3.loc[SDcx.idxmax(), "DF"]+10))*100,0)
    
                    finalLE2=round((round(max(SDc),1)/df4.loc[SDc.idxmax(), "DF"])*100,0)
                    finalLC2=round((round(max(SDc),1)/(df4.loc[SDc.idxmax(), "DF"]+10))*100,0)
    
                    finalLE3=round((round(max(SDc2),1)/df5.loc[SDc2.idxmax(), "DF"])*100,0)
                    finalLC3=round((round(max(SDc2),1)/(df5.loc[SDc2.idxmax(), "DF"]+10))*100,0)
                    
                    with col1: 
                        st.header("Formulation report")
                        st.write("Maximum solubilized drug: "+str(round(max(SDcx),1))+" g/L at "+str(df3.loc[SDcx.idxmax(), "DF"])+" g/L drug feed (LE: "+str(finalLE)+" %, LC: "+str(finalLC)+" %)")
                        max_values = df3.groupby('POL')['SD'].max()
                        max_value = max_values.max()
                        keys_with_max_value = max_values[max_values == max_value].index.tolist()
                        comma_separated_keys = ', '.join(str(key) for key in keys_with_max_value)
                        st.write(comma_separated_keys)
            
                            
                    with col2:
                        im = Draw.MolToImage(Chem.MolFromSmiles(SMI),fitImage=True)
                        st.image(im)
            
                    st.write("Predicted amount of solubilized drug (average by LE and LC models, see below)")
                    fig3=plt.figure(figsize=(10, 6))
                    ax=sns.barplot(x="DF", y="SD", hue="POL", data=df3)
                    plt.xlabel("Drug feed [g/L]")
                    plt.ylabel("Solubilized drug [g/L]")
                    plt.ylim(0, 10)
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    st.pyplot(fig3)
            
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("Amount based on LE models:")                   
                        st.write("Maximum solubilized drug: "+str(round(max(SDc),1))+" g/L at "+str(df4.loc[SDc.idxmax(), "DF"])+" g/L drug feed (LE: "+str(finalLE2)+" %, LC: "+str(finalLC2)+" %)")
                        max_values = df4.groupby('POL')['SD'].max()                    
                        max_value = max_values.max()
                        keys_with_max_value = max_values[max_values == max_value].index.tolist()
                        comma_separated_keys = ', '.join(str(key) for key in keys_with_max_value)
                        st.write(comma_separated_keys)
                        
                        fig1a=plt.figure(figsize=(10, 6))
                        ax = sns.barplot(x="DF", y="SD", hue="POL", data=df4)
                        plt.xlabel("Drug feed [g/L]")
                        plt.ylabel("Solubilized drug [g/L]")
                        plt.ylim(0, 10)
                        ax.get_legend().remove()
                        st.pyplot(fig1a)
            
                        st.write("Calculated from LE predictions:")
                        fig1b=plt.figure(figsize=(10, 6))
                        ax = sns.barplot(x="DF", y="LE", hue="POL", data=df2)
                        plt.xlabel("Drug feed [g/L]")
                        plt.ylabel("Ligand efficiency [%]")
                        plt.ylim(0, 100)
                        ax.get_legend().remove()
                        st.pyplot(fig1b)
                    
                    
                    with col2:
                        st.write("Amount based on LC models:")
                        st.write("Maximum solubilized drug: "+str(round(max(SDc2),1))+" g/L at "+str(df5.loc[SDc2.idxmax(), "DF"])+" g/L drug feed (LE: "+str(finalLE3)+" %, LC: "+str(finalLC3)+" %)")
                        max_values = df5.groupby('POL')['SD'].max()                    
                        max_value = max_values.max()
                        keys_with_max_value = max_values[max_values == max_value].index.tolist()
                        comma_separated_keys = ', '.join(str(key) for key in keys_with_max_value)
                        st.write(comma_separated_keys)
                        fig2a=plt.figure(figsize=(10, 6))
                        ax = sns.barplot(x="DF", y="SD", hue="POL", data=df5)
                        plt.xlabel("Drug feed [g/L]")
                        plt.ylabel("Solubilized drug [g/L]")
                        plt.ylim(0, 10)
                        ax.get_legend().remove()
                        st.pyplot(fig2a)
            
                        st.write("Calculated from LC predictions:")
                        fig2b=plt.figure(figsize=(10, 6))
                        ax = sns.barplot(x="DF", y="LC", hue="POL", data=df2)
                        plt.xlabel("Drug feed [g/L]")
                        plt.ylabel("Loading capacity [%]")
                        plt.ylim(0, 50)
                        ax.get_legend().remove()
                        st.pyplot(fig2b)
            
                    df = pd.read_csv(r'fin_results.csv',index_col=0)
                    df = df.rename(columns={0: "POL", 1:"DF", 2: "LC10", 3: "LC20", 4: "LC30", 5: "LC40", 6: "LE20", 7: "LE40", 8: "LE60", 9: "LE80", 10:"Passed"})
                    df.reset_index(inplace=True)
    
                    st.write("Table of all predictions:")
                    st.dataframe(df.style.applymap(cooling_highlight,subset=["LC10","LC20","LC30","LC40","LE20","LE40","LE60","LE80","Passed"]))     
                    
                    st.write("Amount of passed thresholds:")
                    fig4=plt.figure(figsize=(10, 6))
                    ax=sns.barplot(x="DF", y="Passed", hue="POL", data=df)
                    plt.xlabel("Drug feed [g/L]")
                    plt.ylabel("Thresholds passed")
                    plt.ylim(0, 8)
                    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                    st.pyplot(fig4)

                
                df = pd.read_csv(r'fin_results.csv',index_col=0)
                df = df.rename(columns={0: "POL", 1: "DRUG", 2:"DF", 3: "LC10", 4: "LC20", 5: "LC30", 6: "LC40", 7: "LE20", 8: "LE40", 9: "LE60", 10: "LE80", 11:"Passed"})
                df.reset_index(inplace=True)               
                st.dataframe(df.style.applymap(cooling_highlight,subset=["LC10","LC20","LC30","LC40","LE20","LE40","LE60","LE80","Passed"]))    
    
        for es in ["descriptors.csv","options.csv","create_formulations_temp.R","fin_results.csv","fin_results2.csv","db_test.csv","testformulations.dat","db_formulations.csv"]:
            try:
                os.remove(es)
            except:
                pass

    finally:
        lock.release()

#    except:
#        st.write("Something went wrong. Cannot parse molecules! Please verify your structures.")  

