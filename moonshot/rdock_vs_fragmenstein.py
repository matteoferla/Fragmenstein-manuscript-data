"""
We launch Rdock on Mpro designs and compare then with the results obtained with
Fragmenstein. Fragmenstein results were precomputed and stored in sdf files within the test directory


"""

import os
import shutil
import tempfile

import numpy as np
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem, PandasTools

import chemUtils
from run_rdock import run_rdock, get_rdock_results_sdfname, load_rdock_results
from chemUtils.external.externalToolImporter import ExternalToolImporter
from chemUtils.io import molToFile, molFromFile
from chemUtils.io.io_base import mkdirIfNotExists
from chemUtils.molEdit import remove_atomType
from chemUtils.molEdit.molProtonation import protonate_mols
from chemUtils.scoring.rmsd import calculate_rmsd

import pandas as pd

dirname = os.path.dirname(__file__) #This script parent directory

fragmensteinDir = ExternalToolImporter.get_rootdir("Fragmenstein")
hits_dir = os.path.join(fragmensteinDir, "fragmenstein/mpro/data/hit_mols")
receptor_file = os.path.join(fragmensteinDir, "fragmenstein/mpro/data/template.pdb")

outdir = "/tmp/fragmenstein_dock/mpro/rdock"
mkdirIfNotExists(outdir, useMakedirs=True)

testDir = os.path.join(chemUtils.__path__[0], "..", "tests") #TODO use testtdir instead dirname
fragmenstein_positioned = os.path.join(testDir, "data/positioned.sdf")
fragmenstein_minimised = os.path.join(testDir, "data/minimised.sdf")

n_poses = 15 #Number of poses per compound
timeout = 1000
n_cpus = 1 #Number of compounds to run in parallel
top_K_to_report = 10 #Compute stats with the top_k results

#Load the compounds to be used. Data contains the following columns
#cid,crystal,inspiration_hits
data = pd.read_csv(os.path.join(dirname, 'mproInspirational.csv'), comment="#")

#Load data precalculated with Fragmenstein
fragmenstein_positioned = PandasTools.LoadSDF(fragmenstein_positioned, smilesName='SMILES', molColName='Mol')
fragmenstein_minimised = PandasTools.LoadSDF(fragmenstein_minimised, smilesName='SMILES', molColName='Mol')


#This is when the magic happens. It computes the docking poses using rdock and compare the results with the Fragmenstein
#ones.
#Three different dock runs are conducted, one without constraints, one constrained in which all the pharamacophores
#from the inspirational hits are requested to be satisfied by the solutions and one in which 80% of them are required.
def compute_one_compound(name, cid, fragments):

    resultsDir = os.path.join(outdir, name)
    if os.path.isfile(resultsDir):
        print("%s already computed"%resultsDir)
        df = load_rdock_results(resultsDir)
        return df

    if not os.path.isdir(resultsDir):
        os.mkdir(resultsDir)

    fragments = fragments.split(",")
    fragments_files = [os.path.join(hits_dir, "Mpro-%s.mol"%frag ) for frag in fragments ]
    reference_file = fragments_files[0]
    fragments_mols = [remove_atomType(molFromFile(x)) for x in  fragments_files]

    with tempfile.TemporaryDirectory(suffix=name) as wdir:
        final_outname = get_rdock_results_sdfname(wdir, name=name)[0]



        try:
            mol_for_bond_index = list(fragmenstein_positioned.query("CID == '%s'" % name)['Mol'])[0]
            experimental_mol = molFromFile(os.path.join(hits_dir, "Mpro-"+cid+".mol" ))

            experimental_mol = AllChem.AssignBondOrdersFromTemplate(mol_for_bond_index, experimental_mol)
        except (IndexError, ValueError):
            return

        query_mol = Chem.Mol(experimental_mol)
        query_mol = remove_atomType(query_mol)
        query_mol = protonate_mols([query_mol])[0]
        query_mol = remove_atomType(query_mol)

        ref_mol = molFromFile(reference_file)
        ref_mol = remove_atomType(ref_mol)

        reference_file = os.path.join(wdir, os.path.basename(reference_file))
        molToFile(ref_mol, reference_file)

        df = run_rdock(name=name + "_Freedock", protein=receptor_file, query_mols=query_mol,
                       wdir=os.path.join(wdir, name + "_Freedock"), reference_mol=reference_file, n_poses=n_poses,
                       fragments_for_constraints=None, experimental_mol=experimental_mol, timeout=timeout,
                       verbose=bool(n_cpus == 1))

        df_constraint1 = run_rdock(name + "_Constrained1", receptor_file, query_mol,
                                   os.path.join(wdir, name + "_Constrained1"), reference_file, n_poses=n_poses,
                                   fragments_for_constraints=fragments_mols, fraction_optional_pharm=1.0,
                                   experimental_mol=experimental_mol, timeout=timeout, verbose=bool(n_cpus == 1))
        df_constraint2 = run_rdock(name + "_Constrained2", receptor_file, query_mol,
                                   os.path.join(wdir, name + "_Constrained2"), reference_file, n_poses=n_poses,
                                   fragments_for_constraints=fragments_mols, fraction_optional_pharm=0.8,
                                   experimental_mol=experimental_mol, timeout=timeout, verbose=bool(n_cpus == 1))

        df = pd.concat([df, df_constraint1, df_constraint2], ignore_index=True)

        df = df.reset_index(drop=True)

        fragmenstein_mol = list(fragmenstein_positioned.query("CID == '%s'"%name)['Mol'])[0]
        print(name, cid)

        df["fragmenstein_positioned_RMSD"] = calculate_rmsd(experimental_mol, fragmenstein_mol)

        fragmenstein_mol = list(fragmenstein_minimised.query("CID == '%s'"%name)['Mol'])[0]
        df["fragmenstein_minimised_RMSD"] = calculate_rmsd(experimental_mol, fragmenstein_mol)

        df.sort_values(by="SCORE", inplace=True)
        df["name"] = name

        df = df[["name", 'SMILES', 'Mol', 'SCORE', 'SCORE.norm', "fragmenstein_positioned_RMSD",
                "fragmenstein_minimised_RMSD", "docking_type", "docking_RMSD"]]
        if n_cpus == 1:
            print(df[["name", 'SCORE', 'SCORE.norm', "fragmenstein_positioned_RMSD",
                    "fragmenstein_minimised_RMSD", "docking_type", "docking_RMSD"]].head())
        PandasTools.WriteSDF(df, final_outname,molColName="Mol",properties=list(df.columns))
        shutil.copytree(wdir, resultsDir, dirs_exist_ok=True)
    return df

#Execute the compute_one_compound funtion for all the compouds
list_of_dfs = Parallel(n_jobs=n_cpus, backend="multiprocessing")(delayed(compute_one_compound)(*record) for i, record in data.iterrows())
df = pd.concat(list_of_dfs)

#Save results
PandasTools.WriteSDF(df, os.path.join(outdir, "rdock_vs_fragmenstein_summary.sdf"), molColName="Mol", properties=list(df.columns))
df = df[['name', 'SMILES', 'Mol', 'SCORE', 'SCORE.norm', 'docking_type',
         'docking_RMSD', 'fragmenstein_positioned_RMSD', 'fragmenstein_minimised_RMSD']]
df.to_csv(os.path.join(outdir, "rdock_vs_fragmenstein_summary.csv"), index=False)

#Report statistics
df_ori = df
for to_select_docking in ["pharm4constrdocking", "freedocking"]:
    df = df_ori.query("docking_type == '%s'"%to_select_docking)
    df.reset_index(drop=True, inplace=True)
    def selector_fun(p, k=top_K_to_report):

        p = p.reset_index(drop=True)
        best_score_idxs = np.argsort(p["SCORE"])[:k]
        rmsds = p["docking_RMSD"][best_score_idxs]
        return rmsds.index[rmsds.argmin()]

    selected_idxs = [ part.index[selector_fun(part)] for comp_name, part in df.groupby("name")]

    best_df = df.iloc[selected_idxs,:]
    print("Docking (N=%d): %s"%(len(best_df), to_select_docking))
    for good_thr in [1,2,3]:
        rdock_good = np.sum(best_df['docking_RMSD'] < good_thr)/ best_df.shape[0]
        fragmenstein_good = np.sum(best_df['fragmenstein_minimised_RMSD']< good_thr)/best_df.shape[0]
        print( "RMSD < %f -> Fragmenstein: %f  Rdock: (k = %d) %f" % (good_thr, fragmenstein_good, top_K_to_report, rdock_good))
