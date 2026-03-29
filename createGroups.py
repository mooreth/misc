import math

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

MESH_NSL = 'node.type mesh'
ATOM_NSL= 'node.type atom'
EPS = 1.0e-9                # geometric tolerance

# Use a ray direction that is very unlikely to hit many edges exactly
RAY_DIR = (1.0, 0.3713906763541037, 0.2182178902359924)


# ------------------------------------------------------------
# Small vector utilities
# ------------------------------------------------------------

def v_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def v_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def v_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def v_cross(a, b):
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )

def v_scale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def v_norm(a):
    return math.sqrt(v_dot(a, a))

def v_normalize(a):
    n = v_norm(a)
    if n < EPS:
        return a
    return (a[0]/n, a[1]/n, a[2]/n)


# ------------------------------------------------------------
# Accessors / conversions
# ------------------------------------------------------------

def atom_to_point(atom):
    return (atom.getX().value, atom.getY().value, atom.getZ().value)

def position_to_point(position):
    return (position.x.value, position.y.value, position.z.value)


# ------------------------------------------------------------
# Mesh extraction
# ------------------------------------------------------------

def extract_triangles_from_mesh(mesh):
    """
    Build a flat triangle list [(v0, v1, v2), ...] from all surfaces of a mesh.

    Assumes the positions returned by the surfaces are already in the same
    coordinate frame as the atoms. If your mesh carries a non-identity spatial
    transform, apply that transform to the vertices before classification.
    """
    triangles = []
    transform = mesh.getTransform()

    for surface in mesh.getSurfaces():
        pos = surface.getPositionData()
        idx = surface.getIndexData()

        n_pos = int(len(surface.getPositionData())/3)
        n_tri = int(len(surface.getIndexData())/3)

        verticesPreTransform = [ SBPosition3(SBQuantity.pm(pos[3*i]), SBQuantity.pm(pos[3*i + 1]), SBQuantity.pm(pos[3*i + 2])) for i in range(n_pos) ]
        verticesPostTransform = [ transform.transformPoint(verticesPreTransform[i]) for i in range(n_pos) ]
        vertices = [ position_to_point(verticesPostTransform[i]) for i in range(n_pos) ]

        for t in range(n_tri):
            i0 = int(idx[3*t])
            i1 = int(idx[3*t + 1])
            i2 = int(idx[3*t + 2])
            triangles.append((vertices[i0], vertices[i1], vertices[i2]))

    return triangles


# ------------------------------------------------------------
# Point / triangle predicates
# ------------------------------------------------------------

def point_on_triangle(p, a, b, c, eps=EPS):
    """
    Check whether p lies on triangle abc within tolerance.
    """
    ab = v_sub(b, a)
    ac = v_sub(c, a)
    ap = v_sub(p, a)

    n = v_cross(ab, ac)
    nn = v_dot(n, n)
    if nn < eps:
        return False  # degenerate triangle

    # Distance from point to triangle plane
    dist_num = abs(v_dot(ap, n))
    if dist_num > eps * math.sqrt(nn):
        return False

    # Barycentric test
    d00 = v_dot(ab, ab)
    d01 = v_dot(ab, ac)
    d11 = v_dot(ac, ac)
    d20 = v_dot(ap, ab)
    d21 = v_dot(ap, ac)

    denom = d00 * d11 - d01 * d01
    if abs(denom) < eps:
        return False

    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1.0 - v - w

    return (
        -eps <= u <= 1.0 + eps and
        -eps <= v <= 1.0 + eps and
        -eps <= w <= 1.0 + eps
    )


def ray_intersects_triangle(origin, direction, a, b, c, eps=EPS):
    """
    Moller-Trumbore ray/triangle intersection.
    Returns distance t along the ray if there is a strict hit, else None.

    We exclude hits too close to the origin (t <= eps) to avoid self-counting.
    """
    edge1 = v_sub(b, a)
    edge2 = v_sub(c, a)

    h = v_cross(direction, edge2)
    det = v_dot(edge1, h)

    if abs(det) < eps:
        return None  # parallel or nearly parallel

    inv_det = 1.0 / det
    s = v_sub(origin, a)
    u = inv_det * v_dot(s, h)
    if u < -eps or u > 1.0 + eps:
        return None

    q = v_cross(s, edge1)
    v = inv_det * v_dot(direction, q)
    if v < -eps or u + v > 1.0 + eps:
        return None

    t = inv_det * v_dot(edge2, q)
    if t <= eps:
        return None

    return t


def point_in_watertight_mesh(point, triangles, eps=EPS):
    """
    Parity test:
      - if the point lies on the surface -> treat as inside
      - else count ray/triangle intersections
    """
    # Boundary inclusion
    for a, b, c in triangles:
        if point_on_triangle(point, a, b, c, eps):
            return True

    direction = v_normalize(RAY_DIR)

    # Collect unique hits; this reduces double-counting when the ray passes
    # through shared triangle edges or coplanar patches.
    hits = []

    for a, b, c in triangles:
        t = ray_intersects_triangle(point, direction, a, b, c, eps)
        if t is not None:
            hits.append(t)

    if not hits:
        return False

    hits.sort()

    unique_hits = []
    for t in hits:
        if not unique_hits or abs(t - unique_hits[-1]) > 1.0e-7:
            unique_hits.append(t)

    return (len(unique_hits) % 2) == 1


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

mesh_indexer = SAMSON.getNodes(MESH_NSL)
if len(mesh_indexer) == 0:
    raise RuntimeError(f"No mesh found with NSL: {MESH_NSL}")

mesh = mesh_indexer[0]
print("Mesh:", mesh)

triangles = extract_triangles_from_mesh(mesh)
print("Number of triangles:", len(triangles))

if len(triangles) == 0:
    raise RuntimeError("The mesh contains no triangles.")

atom_indexer = SAMSON.getNodes(ATOM_NSL)

inside_indexer = SBNodeIndexer()
outside_indexer = SBNodeIndexer()

for atom in atom_indexer:
    p = atom_to_point(atom)
    if point_in_watertight_mesh(p, triangles):
        inside_indexer.addNode(atom)
    else:
        outside_indexer.addNode(atom)

print("Atoms inside:", inside_indexer.size)
print("Atoms outside:", outside_indexer.size)

inside_group = SBNodeGroup('inside ' + mesh.name, inside_indexer)
outside_group = SBNodeGroup('outside ' + mesh.name, outside_indexer)

with SAMSON.holding("Create carving groups"):
    inside_group.create()
    SAMSON.getActiveDocument().addChild(inside_group)
    outside_group.create();
    SAMSON.getActiveDocument().addChild(outside_group)
	