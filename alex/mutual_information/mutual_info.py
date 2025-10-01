import mdigest
from mdigest.core.parsetrajectory import *
from mdigest.core.correlation import *
from mdigest.core.dcorrelation import *
from mdigest.core.kscorrelation import *
from mdigest.core.imports import *
import argparse
import warnings
warnings.filterwarnings("ignore")

import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class MutualInfo():
    def __init__(
            self,
            parent : str,
            sr : str,
            lig : str,
            moltype : str,
    ):
        self.parent = Path(parent)
        self.currdir = f'{parent}/{sr}/{moltype}_MI_procedure'
        self.sr = sr
        self.lig = lig
        self.moltype = moltype

        if moltype == "apo":
            self.topology     = f"{parent}/../{moltype}_MDSimulations/50-final-MD-data/md_final_noH2O.gro"
            self.trajectory   = f"{parent}/../{moltype}_MDSimulations/50-final-MD-data/md_cut_aligned_noH2O.xtc"
        elif moltype == "ligand":
            self.topology     = f"{parent}/../{moltype}_MDSimulations/50-final-MD-data/md_{lig}_final_noH2O.gro"
            self.trajectory   = f"{parent}/../{moltype}_MDSimulations/50-final-MD-data/md_{lig}_final_aligned_noH2O.xtc"

        self.savedir = f'{self.currdir}/CACHE'
        mdigest.core.toolkit.folder_exists(self.savedir)

    def processing(self):
        """
        Calculate MI with dynamic correlations, dihedral correlations and electrostatic interactions.
        """
        logging.info(f"Loading {self.moltype} {self.sr} system")
        mds = MDS()
        mds.load_system(self.topology, self.trajectory)
        mds.align_traj(inmem=True, selection='name CA')
        mds.set_selection('protein and name CA', 'protein')
        mds.stride_trajectory(initial=0, final=-1, step=1)


        # correlation from CA displacements 
        logging.info(f"Compute correlation from CA displacements for {self.moltype} {self.sr} system")
        dyncorr = DynCorr(mds)
        dyncorr.parse_dynamics(scale=True, normalize=True, LMI='gaussian', MI='knn_5_2', DCC=True, PCC=True, VERBOSE=True, COV_DISP=True)

        # compute correlation from dihedral fluctuations
        logging.info(f"Compute correlation from dihedral fluctuations for {self.moltype} {self.sr} system")
        dihdyncorr = DihDynCorr(mds)
        dihdyncorr.parse_dih_dynamics(mean_center=True, LMI='gaussian', MI='knn_5_2', DCC=True, PCC=True, COV_DISP=True)

        # Kabsch-Sander analysis for correlation of electrostatic interactions
        logging.info(f"Compute Kabsch-Sander analysis for correlation of electrostatic interactions for {self.moltype} {self.sr} system")        
        KS = KS_Energy(mds)
        # select backbone atoms
        KS.set_selection(['protein and backbone and name N','protein and backbone and name O',
                        'protein and backbone and name C','protein and name H'], system_selstr='protein')
        # set indices according to atom name selection
        KS.set_backbone_dictionary({'N-Backbone':'N',
                                    'O-Backbone':'O',
                                    'C-Backbone':'C',
                                    'CA-Backbone':'CA',
                                    'H-Backbone':'H'})

        KS.KS_pipeline(topology_charges=False, MI="knn_5_2", covariance=True, correction=True)
        KS.compute_EEC(distance_matrix=None, loc_factor=None, don_acc=True)


        # save objects to file for later use
        dyncorr.save_class(file_name_root=f'{self.savedir}/dyncorr')
        KS.save_class(file_name_root=f'{self.savedir}/KS', save_space=True)
        dihdyncorr.save_class(file_name_root=f'{self.savedir}/dih')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate MI using dyncorr, dihcorr, electrostatic interactions.")

    parser.add_argument("-p", "--parent", required=True, type=str, help="Path of parent folder (main SR folder)")
    parser.add_argument("-s", "--sr", required=True, type=str, help="Desired steroid receptor")
    parser.add_argument("-l", "--lig", required=True, type=str, help="Corresponding ligand")
    parser.add_argument("-t","--moltype", required=True, type=str, help="Type of molecule (apo or ligand)")

    args = parser.parse_args()

    parent = Path(args.parent)

    if args.moltype not in ['apo', 'ligand']:
        raise NameError("Invalid moltype, please choose either apo or ligand.")

    mi = MutualInfo(
        parent=str(parent),
        sr = args.sr,
        lig=args.lig,
        moltype=args.moltype
    )

    mi.processing()
