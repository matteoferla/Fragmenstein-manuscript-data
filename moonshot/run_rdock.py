import subprocess
import os
import warnings
from typing import Union, Iterable as typeIterable

import numpy as np
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from rdkit.Chem.Descriptors import MolWt

from chemUtils.external.docking.rdock.pharma4_constraint_extractor import compute_pharm4_constraints
from chemUtils.external.docking.rdock.rdockprm import write_prm_file, convert_ligand_2_file, get_prm_finalname
from chemUtils.external.externalToolImporter import ExternalToolImporter
from chemUtils.io import molToFile, molFromFileOrMol
from chemUtils.molEdit import remove_atomType
from chemUtils.scoring.rmsd import calculate_rmsd
from subprocess import run

_rdock_env = ExternalToolImporter.get_env_vars("rdock")

def _environ_setup(wdir):
    rdock_env = os.environ.copy()
    rdock_env.update(_rdock_env)
    if not os.path.isdir(wdir):
        os.mkdir(wdir)
    rdock_env["RBT_HOME"] = wdir
    return rdock_env

def get_rdock_results_sdfname(wdir, name=""):
    out_file = os.path.join(wdir, "results_"+name)
    final_outname = out_file+"_sorted.sdf"
    return final_outname, out_file

def load_rdock_results(fname):
    df = PandasTools.LoadSDF(fname, smilesName='SMILES', molColName='Mol', strictParsing=False)
    for colname in df.columns:
        if "SCORE" in colname or "RMSD" in colname or "LIGAND_EFFICIENCY" in colname :
            df[colname] = df[colname].astype(float)
    return df

def _prepare_target_rdock(name, reference_mol: Union[str, Chem.Mol], protein: Union[str, Chem.Mol], wdir:str, query_mol=None,
                 fragments_for_constraints=None, timeout=1000, verbose=True,
                          fraction_optional_pharm=1.0, cavity_radius=8):
    """

    :param name: The name of the project
    :param reference_mol: the Chem.Mol or mol file used to compute the pocket
    :param protein: the receptor file or receptor mol
    :param wdir: the directory where results will be saved
    :param query_mol: the query mol if required to compute pharmacophoric constraints
    :param n_poses: the number of poses to produce
    :param fragments_for_constraints: a list of Chem.Mol(s) to extract pharmacophoric restraints for constrained docking
    :param experimental_mol: a  Chem.Mol or mol file molecule to compute rmsd against. Only for testing purposes
    :param timeout: maximum seconds
    :param verbose: whether to print info
    :param fraction_optional_pharm: fraction of pharmacophores to match
    :param cavity_radius: box size radius
    :return:
    """
    wdir = os.path.abspath(wdir)
    rdock_env =  _environ_setup(wdir)

    ref_mol = molFromFileOrMol(reference_mol)
    print(Chem.MolToSmiles(ref_mol))

    ref_mol = remove_atomType(ref_mol)
    reference_file= convert_ligand_2_file(ref_mol, wdir=wdir)

    pharma_constraints_required, pharma_constraints_optional = None, None
    if fragments_for_constraints is not None:
        assert query_mol is not None, "Error, if setting pharmacophoric constraints, you need to provide a query molecule"
        query_mol = molFromFileOrMol(query_mol)
        if query_mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(query_mol)
        fragments_for_constraints = [ molFromFileOrMol(molOrFname) for molOrFname in fragments_for_constraints ]

        (pharma_constraints_required, pharma_constraints_optional,
         n_query_pharm) = compute_pharm4_constraints(fragments_for_constraints, query_mol=query_mol, use_mcs=True)
        n_optional_pharm = len(pharma_constraints_optional)
        if n_query_pharm is not None:
            n_optional_pharm = min(n_optional_pharm, n_query_pharm)

        n_optional_pharm = max(1, int(n_optional_pharm*fraction_optional_pharm))
        phar4_weight = 1 + abs(1-fraction_optional_pharm)
    else:
        n_optional_pharm = 1
        phar4_weight = 1
    params_file = write_prm_file(name=name, wdir=wdir, title=name, receptor_file=protein,
                                 cavity_method="REFERENCE_LIGAND",
                                 reference_ligand_file=reference_file, cavity_radius=cavity_radius,
                                 pharma_constraints_required=pharma_constraints_required,
                                 pharma_constraints_optional = pharma_constraints_optional,
                                 num_optional_pharma_to_match=n_optional_pharm, phar4_weight=phar4_weight)
    params_file = os.path.abspath(params_file)
    cmd = [os.path.join(_rdock_env["RDOCK_BINDIR"], 'rbcavity'), '-was', '-d', '-r', params_file]
    cmd = list(map(str, cmd))
    if verbose:
        print( " ".join([x+"="+y for x,y in _rdock_env.items()]), " ".join(cmd))

    out = run(cmd, cwd=wdir,  timeout=timeout, check=True, capture_output= True, env=rdock_env)
    assert "DOCKING SITE" in out.stdout.decode("utf-8"), "%s: Error, cavity preparation failed for %s"%(  out.stdout.decode("utf-8"), name)
    if verbose:
        print( out.stdout.decode("utf-8"))

    return params_file

def _launch_rdock(name, query_mol, wdir, params_file, n_poses=10, experimental_mol=None, timeout=1000, verbose=False):
    # print(Chem.MolToSmiles(query_mol))
    wdir = os.path.abspath(wdir)
    final_outname, out_file = get_rdock_results_sdfname(wdir, name=name)
    if os.path.isfile(final_outname):
        print("%s already computed"%final_outname)
        df = load_rdock_results(final_outname)
        return df

    rdock_env =  _environ_setup(wdir)
    query_mol_fname = os.path.join(wdir, name + ".sdf")

    if not os.path.isfile(query_mol_fname):
        query_mol = molFromFileOrMol(query_mol)
        if query_mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(query_mol)

        molToFile(query_mol, query_mol_fname)

    cmd = [os.path.join(_rdock_env["RDOCK_BINDIR"], 'rbdock'), '-r', os.path.abspath(params_file), "-p", "dock.prm",
         "-n", n_poses, "-i", os.path.abspath(query_mol_fname), "-o", os.path.abspath(out_file)]
    cmd = list(map(str, cmd))
    if verbose:
        print( " ".join([x+"="+y for x,y in _rdock_env.items()]) + " " + " ".join(cmd))
    try:
        out2 = run(cmd, cwd=wdir, timeout=timeout, check=True, capture_output=True, env = rdock_env)
    except  subprocess.CalledProcessError:
        print("Unexpected error launching %s"%(" ".join(cmd)))
        raise
    if verbose:
        print( out2.stdout.decode("utf-8"))
    if  "_ERROR" in out2.stdout.decode("utf-8"):
        lines =  out2.stdout.decode("utf-8").split("\n")
        for line in lines:
            if "_ERROR"  in line:
                assert line.startswith("**WARNING**"), "%s: Error running docking %s"%(
                                                                                out2.stdout.decode("utf-8"), name)
    sdf_fname = out_file+".sd"
    assert os.path.isfile(sdf_fname), "%s: Error running docking %s" % (
            out2.stdout.decode("utf-8"), name)
    # print(out.stdout.decode("utf-8"))

    df = load_rdock_results(sdf_fname)

    useful_colnames = ["name", 'SMILES', 'Mol', 'SCORE', 'SCORE.norm']

    if experimental_mol is not None:
        experimental_mol = molFromFileOrMol(experimental_mol)
        df["docking_RMSD"] = df["Mol"].map( lambda mol: calculate_rmsd(experimental_mol, mol) )
    else:
        df["docking_RMSD"] = np.nan

    with open(params_file) as f:
        if "PHARMACOPHORIC RESTRAINTS" in f.read():
            df["docking_type"] = "pharm4constrdocking"
        else:
            df["docking_type"] = "freedocking"

    useful_colnames += ["docking_type", "docking_RMSD"]


    df.reset_index(drop=True, inplace=True)

    df["name"] = name

    useful_colnames += ["LIGAND_EFFICIENCY", "LIGAND_EFFICIENCY.norm"]

    if len(df) == 0:
        for useful_col in useful_colnames:
            df[useful_col] = np.nan
    else:
        df["LIGAND_EFFICIENCY"] = df.apply(lambda row: float(row["SCORE"])/MolWt(row["Mol"]), axis=1)
        df["LIGAND_EFFICIENCY.norm"] = df.apply(lambda row: float(row["SCORE.norm"])/MolWt(row["Mol"]), axis=1)
        df = df[useful_colnames]

    PandasTools.WriteSDF(df, final_outname,molColName="Mol",properties=list(df.columns))
    df.sort_values(by="SCORE", inplace=True)
    return df

def run_rdock(name, protein: Union[str, Chem.Mol], query_mols: Union[str, Chem.Mol, typeIterable[Union[str, Chem.Mol]]],
              wdir: str, reference_mol: Union[str, Chem.Mol], n_poses=10, n_cpus: int = 1,
              fragments_for_constraints=None, fraction_optional_pharm=1.0, experimental_mol=None, timeout=1000,
              verbose=False):
    """
    Executes a docking run
    :param name: The name of the docking experiment
    :param protein: the protein (receptor) filename or Chem.Mol
    :param query_mols: The query molecule(s) or filename(s) to dock
    :param wdir: Directory were files are going to be saved
    :param reference_mol: A reference mol or filaname to compute the binding site. E.g. a different bindier with structure
    :param n_poses: Number of poses to calculate
    :param n_cpus: Number of cpus to use in parallel. Paralleization happens at the query_mol vlevel
    :param fragments_for_constraints: A set of fragments to extracted pharmacophoric constraints before docking
    :param fraction_optional_pharm: Fraction (0.5-1) of pharmacophores to be required to be satisfied.
    :param experimental_mol: An experimental ground truth pose of the query_molecule. Only for benchmarking purposes
    :param timeout: Timeout in seconds to interrumpt execution
    :param verbose: Print a lot of logs to stdin
    :return:
    """
    if not isinstance(query_mols, (list, tuple)):
        query_mols = [query_mols]
        isIterable = False
    else:
        query_mols = query_mols
        isIterable = True

    results = []
    if fragments_for_constraints:
        for i, query_mols in enumerate(query_mols):
            params_file = get_prm_finalname(name+"_%d"%i, wdir)
            if not os.path.exists(params_file):
                params_file = _prepare_target_rdock(name=name + "_%d" % i, reference_mol=reference_mol, protein=protein,
                                                    wdir=wdir, query_mol=query_mols,
                                                    fragments_for_constraints=fragments_for_constraints,
                                                    timeout=timeout, verbose=verbose,
                                                    fraction_optional_pharm=fraction_optional_pharm)
            df = _launch_rdock(name=name, query_mol=query_mols, wdir=wdir, params_file=params_file, n_poses=n_poses,
                               experimental_mol=experimental_mol, timeout=timeout, verbose=verbose)
            results.append(df)
    else:
        params_file = get_prm_finalname(name, wdir)
        if not os.path.exists(params_file):
            params_file = _prepare_target_rdock(name=name, reference_mol=reference_mol, protein=protein, wdir=wdir,
                                                query_mol=query_mols,
                                                fragments_for_constraints=fragments_for_constraints, timeout=timeout,
                                                verbose=verbose)
        def fun(i, query_mol):
            return _launch_rdock(name=name + "_%d" % i, query_mol=query_mol, wdir=wdir, params_file=params_file,
                                 n_poses=n_poses, experimental_mol=experimental_mol, timeout=timeout, verbose=verbose)

        results = Parallel(n_jobs=n_cpus)(delayed(fun)( i, query_mol) for i, query_mol in enumerate(query_mols) )

    if not isIterable:
        return results[0]
    else:
        return results
