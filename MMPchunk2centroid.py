import re

def print_chunk_centroids(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    current_chunk = None
    atoms = []

    for line in lines:

        # --------------------------------------------------
        # Look for the beginning of a chunk
        # mol (ChunkName) def
        # --------------------------------------------------

        mol_match = re.search(r'mol\s+\((.*?)\)\s+def', line)

        if mol_match:

            # Print the previous chunk before starting the new one
            if current_chunk is not None and atoms:
                print_chunk_centroid(current_chunk, atoms)

            current_chunk = mol_match.group(1)
            atoms = []

            continue

        # --------------------------------------------------
        # Look for atom definitions
        # --------------------------------------------------

        if line.find("atom ") > -1 and line.find("def") > -1:

            # Remove parentheses and commas
            clean_line = line.replace('(', '').replace(')', '')
            clean_line = clean_line.replace(',', '')

            mmp_atoms = clean_line.split()

            # MMP format:
            # atom INDEX PROTON_NUMBER X Y Z ...

            try:
                x = float(mmp_atoms[3])
                y = float(mmp_atoms[4])
                z = float(mmp_atoms[5])

                atoms.append((x, y, z))

            except (IndexError, ValueError):
                print("Could not parse atom:", line.strip())

    # ------------------------------------------------------
    # Print the final chunk
    # ------------------------------------------------------

    if current_chunk is not None and atoms:
        print_chunk_centroid(current_chunk, atoms)


def print_chunk_centroid(chunk_name, atoms):

    num_atoms = len(atoms)

    x_centroid = sum(atom[0] for atom in atoms) / num_atoms /10000
    y_centroid = sum(atom[1] for atom in atoms) / num_atoms /10000
    z_centroid = sum(atom[2] for atom in atoms) / num_atoms /10000

    print(
        f"{chunk_name}: "
        f"{num_atoms} atoms, "
        f"centroid = "
        f"({x_centroid:.6f}, "
        f"{y_centroid:.6f}, "
        f"{z_centroid:.6f})"
    )


# ----------------------------------------------------------
# Example
# ----------------------------------------------------------

filename = "D:/LAMMPS projects/Clutch/clutch_assembly.mmp"

print_chunk_centroids(filename)