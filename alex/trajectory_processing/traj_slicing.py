import mdtraj as md
import numpy as np

#------------------------------------------------------------------------------------------------#

def folded_unfolded_division(desired_traj, rmsd_threshold):
    """
    Function used to separate the folded and the unfolded part from a given trajectory according to RMSD.
    The function needs as input the trajectory as `.xtc` loaded with MDtraj and the RMSD threshold.
    Returns a tuple with the folded and unfolded trajectories.
    """
    # computing RMSD 
    rmsd_partialtraj = md.rmsd(desired_traj, first_frame_ref)
    rmsd_angstrom = rmsd_partialtraj*10

    # getting indices for folded/unfolded part
    FOLDED_part_ndx = np.where(rmsd_angstrom <= rmsd_threshold)[0]
    UNFOLDED_part_ndx = np.where(rmsd_angstrom > rmsd_threshold)[0]

    FOLDED_part_traj = partial_traj[FOLDED_part_ndx]
    UNFOLDED_part_traj = partial_traj[UNFOLDED_part_ndx]

    return FOLDED_part_traj, UNFOLDED_part_traj
    
#------------------------------------------------------------------------------------------------#

def traj_pruning(sliced_part_traj, desired_N_frames):
    """
    Function used to reduce the frame number from a given trajectory to be suitable for excogito optimisation.
    It takes as input any trajectory loaded with MDtraj in the format `.xtc`, and the number of frames.
    The function returns the pruned trajectory, removing frames with a stride equal to the ratio between the total number of frames and the desired number of frames.
    If the trajectory length is already smaller than the desired number of frames, then the same unpruned trajectory is returned.
    """
    tot_frames = len(sliced_part_traj)
    print(f"Tot frames in the current trajectory: {tot_frames}")
    if tot_frames > desired_N_frames:
        stride = tot_frames / desired_N_frames
        print(f"Stride: {stride}")
        frame_ndx = np.arange(0, tot_frames, stride)
        frame_ndx=frame_ndx.astype(int)
        reduced_traj = sliced_part_traj[frame_ndx]
        print(f"Length of new trajectory: {len(reduced_traj)}")
        return reduced_traj
    else:
        print(f"Your trajectory already satisfies the number of frames: {tot_frames} < 10000")
        return sliced_part_traj

#------------------------------------------------------------------------------------------------#

