import crystals
import numpy as np
import os

from sys import argv


# FUNCTION DEFINITIONS

def all_atoms(crystal):
    return [ a for a in sorted(crystal.supercell(3, 3, 3)) ]


def int_atoms(crystal, which='Al'):
    is_in_middle = lambda x: all(1 <= i <= 2 for i in x)

    return [ a for a in sorted(crystal.supercell(3, 3, 3)) if a.symbol == which and is_in_middle(a.coords_fractional) ]


def lowest_distances(atom, atoms_other):
    names_other = [ a.symbol for a in atoms_other ]
    distances = [ crystals.distance_cartesian(atom, o) for o in atoms_other if crystals.distance_cartesian(atom, o) != 0. ]
    radii = [ atom.coords_cartesian - o.coords_cartesian for o in atoms_other if crystals.distance_cartesian(atom, o) != 0. ]

    return sorted(zip(names_other, distances, radii), key=lambda x: x[1])[:]


def angle(r1, r2):
    return np.arccos(np.dot(r1, r2)/np.sqrt(np.dot(r1, r1) * np.dot(r2, r2))) * 180/np.pi


def in_range(value, midpoint, width):
    return midpoint - width < value < midpoint + width


# CONSTANTS
ATOM_OF_INTEREST = 'Al'
DISTANCE_TOL_C4 = 1.2
DISTANCE_TOL_C6 = 1.2
ANGLE_TOL_C4 = 20
ANGLE_TOL_C6_MIN = 20
ANGLE_TOL_C6_MAJ = 25
ANGLE_VAL_C4 = 109.5
ANGLE_VAL_C6_MIN = 90
ANGLE_VAL_C6_MAJ = 180

# MAIN BODY
if __name__ == '__main__':
    try:
        filename = argv[1]
    except:
        print('Expected filename argument.')
        exit(1)

    try:
        c = crystals.Crystal.from_cif(filename)
    except Exception as e:
        print(e)
        print('File not found.')
        exit(1)

    prefix = filename.split('/')[-1].split('.')[0]

    os.makedirs(prefix, exist_ok=True)

    aa = all_atoms(c)
    ias = int_atoms(c, which=ATOM_OF_INTEREST)
    coordinations = [0] * len(ias)

    content = 'element,x,y,z\n' + \
        '\n'.join(f'{a.symbol},{a.coords_cartesian[0]},{a.coords_cartesian[1]},{a.coords_cartesian[2]}' for a in aa)

    with open(f'{prefix}/{prefix}_all_atoms_supercell.csv', 'w') as f:
        f.write(content)

    for n, ia in enumerate(ias):
        LEN = 8
        sa = lowest_distances(ia, aa)[:LEN]  # surrounding LEN atoms

        # write atom number, coordinates in debug file
        content = 'index,element,distance,dx,dy,dz\n' + \
            '\n'.join(f'{i},{name},{d},{r[0]},{r[1]},{r[2]}' for i, (name, d, r) in enumerate(sa))

        with open(f'{prefix}/{prefix}_{ATOM_OF_INTEREST}{str(n).zfill(2)}_closest{LEN}.csv', 'w') as f:
            f.write(content)

        # write angles in debug file:
        index_pairs = [ (i, j) for i in range(LEN) for j in range(i+1, LEN) ]
        angles_pairwise = { (i, j): angle(sa[i][2], sa[j][2]) for i, j in index_pairs }

        content = ' \t' + '\t'.join(str(i) for i in range(LEN)) + '\n' + \
            '\n'.join((f'{i}\t' + ' \t'*i + 'x\t' + '\t'.join(f'{angles_pairwise[(i, j)]:.2f}' for j in range(i+1, LEN))) for i in range(LEN));

        with open(f'{prefix}/{prefix}_{ATOM_OF_INTEREST}{str(n).zfill(2)}_angles.csv', 'w') as f:
            f.write(content)

        # octohedral coordination check
        angles_oct = sorted([ a for (i, j), a in angles_pairwise.items() if i < 6 and j < 6 ])

        distance_cond_oct = all(a[1] < sa[0][1] * DISTANCE_TOL_C6 for a in sa[:6])
        angle_cond_oct = all(in_range(a, ANGLE_VAL_C6_MIN, ANGLE_TOL_C6_MIN) for a in angles_oct[:-3]) and \
            all(in_range(a, ANGLE_VAL_C6_MAJ, ANGLE_TOL_C6_MAJ) for a in angles_oct[-3:])

        if distance_cond_oct and angle_cond_oct:
            coordinations[n] = 6
            print(f'{ATOM_OF_INTEREST}{str(n).zfill(2)} - octohedral coordination')
            continue

        # tetrahedral coordination check
        angles_tet = sorted([ a for (i, j), a in angles_pairwise.items() if i < 4 and j < 4 ])

        distance_cond_tet = all(a[1] < sa[0][1] * DISTANCE_TOL_C4 for a in sa[:4])
        angle_cond_tet = all(in_range(a, ANGLE_VAL_C4, ANGLE_TOL_C4) for a in angles_tet)

        if distance_cond_tet and angle_cond_tet:
            coordinations[n] = 4
            print(f'{ATOM_OF_INTEREST}{str(n).zfill(2)} - tetrahedral coordination')
            continue

        coordinations[n] = 0
        print(f'{ATOM_OF_INTEREST}{str(n).zfill(2)} - no coordination found')

    content = 'index,coordination\n' + '\n'.join(f'{i},{n}' for i, n in enumerate(coordinations))
    with open(f'{prefix}_coordination.csv', 'w') as f:
        f.write(content)
