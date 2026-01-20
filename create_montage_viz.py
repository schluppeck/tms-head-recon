# create_montage_viz.py
# ds 2026-01-14, create a simple viz of electrode montage on head surface
#
# with quite a bit of help from Claude in VS code... getting more and more useful!

# on my machine, I use conda env for MNE
# If using VS Code, set interpreter to:

# /Applications/MNE-Python/1.11.0_0/.mne-python/bin/python
import os
import sys
import numpy as np
import mne
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from pathlib import Path


# check command line args
if len(sys.argv) == 2:
    prefer = sys.argv[1]
    if prefer not in ['up', 'down']:
        print("Usage: python create_montage_viz.py [up|down]")
        sys.exit(1)
else:
    prefer = 'up'  # default preference


def key_and_exit():
    """Wait for user to press enter key before exiting

    Simple function to wait for user keypress before exiting
    """
    print("Press the enter key to close window...")
    input()
    sys.exit(0)


def build_adjacency(vertices, faces):
    """Generate adjacency matrix from mesh vertices and faces

       Diijkstra algorithm needs adjacency matrix to compute geodesic distances 
    """
    n = len(vertices)
    edges = []
    weights = []

    for tri in faces:
        for i in range(3):
            v1, v2 = tri[i], tri[(i+1) % 3]
            dist = np.linalg.norm(vertices[v1] - vertices[v2])
            edges.append([v1, v2])
            weights.append(dist)

    edges = np.array(edges)
    adj = csr_matrix((weights, (edges[:, 0], edges[:, 1])), shape=(n, n))
    return adj + adj.T  # make symmetric


def reconstruct_path(pred_matrix, source_idx, dest_idx, source_row=0):
    """Reconstruct shortest path from predecessor matrix

    Args:
        pred_matrix: Predecessor matrix from dijkstra
        source_idx: Index of source vertex
        dest_idx: Index of destination vertex
        source_row: Row in pred_matrix corresponding to the source (0 for nasion, 1 for inion)

    Returns:
        List of vertex indices from source to destination
    """
    path = [dest_idx]
    current = dest_idx

    while current != source_idx:
        current = int(pred_matrix[source_row, current])
        if current == -9999:  # dijkstra uses this as invalid marker
            return None
        path.append(current)

    return path[::-1]  # reverse to get source->dest order


def build_sagittal_constrained_adjacency(vertices, faces, x_tolerance=5, z_prefer='down'):
    """Generate adjacency matrix with penalty for deviating from sagittal plane

    This constrains paths to stay near a constant x-value by adding high costs
    to edges that cross far from the median x-coordinate.

    Args:
        vertices: Mesh vertices in mm
        faces: Mesh face indices
        x_tolerance: How much x-deviation (mm) is allowed before high penalty

    Returns:
        Constrained adjacency matrix
    """
    n = len(vertices)
    edges = []
    weights = []

    # Target x-range (median of all x values)
    x_median = np.median(vertices[:, 0])

    for tri in faces:
        for i in range(3):
            v1, v2 = tri[i], tri[(i+1) % 3]
            dist = np.linalg.norm(vertices[v1] - vertices[v2])

            # Add penalty for x-deviation from median
            x_dev = (abs(vertices[v1, 0] - x_median) +
                     abs(vertices[v2, 0] - x_median)) / 2

            # high penalty for deviation
            penalty = max(0, (x_dev - x_tolerance) * 10)

            # Add penalty if path has z values that are either negative or positive
            z_pos = vertices[v2, 2] + vertices[v1, 2]

            if z_prefer == 'up' and z_pos < 0:
                penalty += abs(z_pos) * 10

            # reality check: prefer downward paths?
            if z_prefer == 'down' and z_pos >= 0:
                penalty += abs(z_pos) * 10

            total_weight = dist + penalty
            edges.append([v1, v2])
            weights.append(total_weight)

    edges = np.array(edges)
    adj = csr_matrix((weights, (edges[:, 0], edges[:, 1])), shape=(n, n))
    return adj + adj.T  # make symmetric


SUBJECTS_DIR = os.getenv("SUBJECTS_DIR")
sub = "subject-A"

brain = mne.viz.Brain(
    subject=sub,
    hemi="both",
    surf="pial",
    subjects_dir=SUBJECTS_DIR,
    background="white",
    size=(800, 600),
)

# FEF-ish
labels = (
    "ctx-lh-precentral",
    "ctx-lh-superiorfrontal",
    "ctx-rh-precentral",
    "ctx-rh-superiorfrontal",
)

brain.add_volume_labels(aseg="aparc+aseg", labels=labels)


# these are 3d points (taken in mm from freeview window)
# kind of by hand
inion_coords = np.array([-14, -95, 12])  # in millimeters
nasion_coords = np.array([3, 92, -3.5])  # in millimeters


skin_surface_file = Path(SUBJECTS_DIR) / sub / "bem" / "subject-A-head.fif"
print("Loading skin surface from:", skin_surface_file)
# Load BEM surface (skin surface is usually the outermost)
bem_surfaces = mne.read_bem_surfaces(skin_surface_file)

# skin is typically index 2 (0=inner skull, 1=outer skull, 2=skin)
skin_surface = bem_surfaces[0]

# # Get the vertices and faces
vertices = skin_surface['rr']  # vertices in meters
faces = skin_surface['tris']  # triangle indices

# # Manually identify nasion and inion indices on the surface
nasion_idx = 4826  # for subject-A
inion_idx = 7818

# reality check does the IDX lead to the right coords?
nasion_coords_from_mesh = vertices[nasion_idx] * 1000  # to mm
inion_coords_from_mesh = vertices[inion_idx] * 1000  # to mm

# you can check if they are close
np.isclose(nasion_coords, nasion_coords_from_mesh, rtol=1e-1)

# plot them to see...
points = np.array([nasion_coords, inion_coords])  # nasion and inion
# add foci, red, scale factor 1mm,
brain.add_foci(points,
               color='blue',
               scale_factor=1,
               hemi="vol")

# project onto the pial surface
brain.add_foci(points,
               color='red',
               scale_factor=1,
               map_surface="pial",
               hemi="lh")

brain.add_head(dense=False)

brain.show_view(azimuth=90, elevation=0, distance=800)  # 500mm

# get adjacency matrix

# make sure vertices are in mm
vertices_mm = vertices * 1000

# Use sagittal-constrained adjacency for path along constant x
adj_matrix = build_sagittal_constrained_adjacency(
    vertices_mm, faces, x_tolerance=5, z_prefer=prefer)  # pass up/down preference

# compute geodesic with sagittal constraint
distances, pd = dijkstra(adj_matrix, indices=[
                         nasion_idx, inion_idx], return_predecessors=True)


# Find path from nasion to inion
path_nasion_to_inion = reconstruct_path(
    pd, nasion_idx, inion_idx, source_row=0)
print(f"Path from nasion to inion: {path_nasion_to_inion}")
print(f"Path length (num vertices): {len(path_nasion_to_inion)}")
print(f"Geodesic distance: {distances[0, inion_idx]:.2f} mm")

# Find midpoint along the path
# initalize approximate method (not quite correct)
# half the way by index... proper way: by distnace... but more code!
midpoint_idx = path_nasion_to_inion[len(path_nasion_to_inion) // 2]
midpoint_coords = vertices_mm[midpoint_idx]

print(f"\nMidpoint vertex index: {midpoint_idx}")
print(f"Midpoint coordinates (mm): {midpoint_coords}")

# plot them...
path_points = vertices[path_nasion_to_inion] * 1000  # to mm
brain.add_foci(path_points,
               color='green',
               scale_factor=0.2,
               hemi="vol")

# Plot midpoint in yellow
brain.add_foci(midpoint_coords,
               color='yellow',
               scale_factor=1.5,
               hemi="vol")

# Plot midpoint PROJECTED onto the pial surface!


brain.add_foci([midpoint_coords],
               color='yellow',
               scale_factor=1.5,
               map_surface="pial",
               hemi='lh')

# shortest path is green dots ... but need to enforce this to be in the sagittal plan

key_and_exit()
