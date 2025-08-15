import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_shape(ax, verts, faces, title):
    ax.add_collection3d(Poly3DCollection([verts[face] for face in faces], alpha=0.6, edgecolor='k'))
    ax.scatter3D(*verts.T, color='r')
    ax.set_title(title)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)

def make_reference_prism():
    # Two triangle faces (top & bottom)
    h = np.sqrt(3.) / 4
    top = np.array([[0, 1, h], [0, 0, h], [1, 0, h]])
    bottom = top.copy()
    bottom[:, 2] = -1.
    verts = np.vstack([top, bottom])
    faces = [
        [0, 1, 2],  # top
        [3, 4, 5],  # bottom
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [2, 0, 3, 5]
    ]
    return verts, faces

def make_target_prism():
    # Two triangle faces (top & bottom)
    h = np.sqrt(3) / 2
    top = np.array([[0, 1, h], [-np.sqrt(3)/2, -0.5, h], [np.sqrt(3)/2, -0.5, h]])
    bottom = top.copy()
    bottom[:, 2] = 0
    verts = np.vstack([top, bottom])
    faces = [
        [0, 1, 2],  # top
        [3, 4, 5],  # bottom
        [0, 1, 4, 3],
        [1, 2, 5, 4],
        [2, 0, 3, 5]
    ]
    return verts, faces

def make_tetrahedron():
    # Vertices of a regular tetrahedron centered at origin
    verts = np.array([
        [1, 0, -1/np.sqrt(2)],
        [-0.5, np.sqrt(3)/2, -1/np.sqrt(2)],
        [-0.5, -np.sqrt(3)/2, -1/np.sqrt(2)],
        [0, 0, np.sqrt(2)]
    ])
    faces = [
        [0, 1, 2],
        [0, 1, 3],
        [1, 2, 3],
        [2, 0, 3]
    ]
    return verts, faces

def make_square_pyramid():
    # Square base pyramid
    verts = np.array([
        [-1, -1, 0],
        [1, -1, 0],
        [1, 1, 0],
        [-1, 1, 0],
        [0, 0, 1.5]
    ])
    faces = [
        [0, 1, 2, 3],  # base
        [0, 1, 4],
        [1, 2, 4],
        [2, 3, 4],
        [3, 0, 4]
    ]
    return verts, faces


# reference Prism
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(111, projection='3d')
verts, faces = make_reference_prism()
plot_shape(ax1, verts, faces, "Reference Prism")

# target Prism
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(111, projection='3d')
verts, faces = make_target_prism()
plot_shape(ax1, verts, faces, "Target Prism")

## Tetrahedron
#fig = plt.figure(figsize=(7, 7))
#ax2 = fig.add_subplot(111, projection='3d')
#verts, faces = make_tetrahedron()
#plot_shape(ax2, verts, faces, "Tetrahedron")
#
## Pyramid
#fig = plt.figure(figsize=(7, 7))
#ax3 = fig.add_subplot(111, projection='3d')
#verts, faces = make_square_pyramid()
#plot_shape(ax3, verts, faces, "Square Pyramid")

plt.tight_layout()
plt.show()

