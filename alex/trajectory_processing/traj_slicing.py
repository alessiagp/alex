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
    Reduce the number of frames in an MDtraj trajectory to exactly desired_N_frames.

    Parameters
    ----------
    sliced_part_traj : mdtraj.Trajectory
        The input trajectory.
    desired_N_frames : int
        Target number of frames after pruning.

    Returns
    -------
    reduced_traj : mdtraj.Trajectory
        The pruned trajectory. If the input has fewer frames than desired_N_frames,
        the original trajectory is returned unchanged.
    """
    tot_frames = len(sliced_part_traj)
    print(f"Total frames in the current trajectory: {tot_frames}")

    if tot_frames > desired_N_frames:
        # Create evenly spaced indices, always exactly desired_N_frames
        frame_ndx = np.round(
            np.linspace(0, tot_frames - 1, desired_N_frames)
        ).astype(int)

        # Safety check: remove duplicates if rounding causes overlaps
        frame_ndx = np.unique(frame_ndx)

        # If we lost some frames due to duplicates, add extras from the end
        while len(frame_ndx) < desired_N_frames:
            frame_ndx = np.append(frame_ndx, tot_frames - 1)
            frame_ndx = np.unique(frame_ndx)

        reduced_traj = sliced_part_traj[frame_ndx]
        print(f"Reduced trajectory length: {len(reduced_traj)}")
        return reduced_traj
    else:
        print(f"Trajectory already has {tot_frames} frames, which is ≤ desired {desired_N_frames}.")
        return sliced_part_traj

#------------------------------------------------------------------------------------------------#

