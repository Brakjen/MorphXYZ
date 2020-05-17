import numpy as np
import os
import json


def get_atomic_mass(sym):
    """
    Return the atomic mass of passed atomic symbol. Myst be IUPAC symbol (case sensitive).
    :param sym:
    :return: float
    """
    with open("atomic_masses.json") as f:
        m = json.load(f)
    return float(m[sym])


def get_coords(xyzfile):
    with open(xyzfile) as f:
        tmp = f.readlines()[2:]
        tmp = [line.split() for line in tmp]
        labels = [line[0] for line in tmp]
        tmp = [[float(c) for c in line[1:]] for line in tmp]
        return np.asarray(tmp), np.asarray(labels)


def com(P, labels):
    masses = np.asarray([get_atomic_mass(atom) for atom in labels])
    return np.einsum("ki, k -> i", P, masses) / sum(masses)


def center_centroid(A):
    """
    Return array whose centroid are at the origin.
    :param A:
    :param B:
    :return: np.ndarray
    """
    return A - centroid(A)


def center_com(A, l):
    return A - com(A, l)


def centroid(A):
    """
    Return the centroid of the passed array A.
    :param A: Array for which centroid should be computed
    :return: np.ndarray
    """
    return np.average(A, axis=0)


def covariance(A, B):
    """
    Return the covariance matrix of A and B.
    :param A:
    :param B:
    :return: np.Array
    """
    return np.einsum("ki, kj -> ij", A, B)


def write_xyz(coords, labels, fname):
    """
    Write XYZ coordinates to file
    :param coords:
    :param labels:
    :param fname: output filename
    :return:
    """
    with open(fname, "w") as f:
        f.write(f"{coords.shape[0]}\n")
        f.write("\n")
        for atom, label in zip(coords, labels):
            f.write(f"{label} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n")


def optimal_rotation(A, B, verbose=False, method="centroid", la=None, lb=None):
    """
    Return the optimal rotation matrix.
    :param C: Covariance matrix
    :return:
    """
    # Align molecules
    if method == "centroid":
        A, B = center_centroid(A), center_centroid(B)
        if verbose:
            print("Origin start:", centroid(A))
            print("Origin stop :", centroid(B))
    elif method == "com":
        A, B = center_com(A, la), center_com(B, lb)
        if verbose:
            print("Origin start", com(A, la))
            print("Origin stop ", com(B, lb))

    # Covariance matrix
    C = covariance(A, B)
    if verbose:
        print("Covariance Matrix:")
        print(C)

    # Singular Value Decomposition
    U, s, V = np.linalg.svd(C)

    # Compute determinant
    d = np.linalg.det(np.einsum("ik, kj -> ij", U, V))

    # Correction matrix
    K = np.identity(3) * np.array([1, 1, d])
    if verbose:
        print("Correction Matrix:")
        print(K)
    return np.einsum("ik, kl, lj -> ij", U, K, V)


def align(P, Q, verbose=False, method="centroid", la=None, lb=None):
    """
    Align mol1 onto mol2.
    :param mol1: XYZ file. Molecule to be aligned.
    :param mol2: XYZ file. Reference molecule.
    :param verbose: print matrices
    :return: np.Array
    """
    R = optimal_rotation(P, Q, verbose=verbose, method=method, la=la, lb=lb)
    if verbose:
        print("Optimal Rotation Matrix:")
        print(R)

    new = np.einsum("ik, kj -> ij", P, R)
    return new - centroid(new) + centroid(Q)


def interpolate(A, B, nsteps=100, verbose=False):
    """"""
    if verbose:
        print("")
        print("Number of steps: ", nsteps)
        print("Aligning...")

    A = align(A, B)
    interp_space = np.linspace(0, 1, nsteps)

    if verbose:
        print("Interpolating...")
    forward = np.array([np.zeros((A.shape[0], 3)) for i in range(nsteps)])
    backward = np.array([np.zeros((A.shape[0], 3)) for i in range(nsteps)])

    for i, scale in enumerate(interp_space):
        forward[i] = A + scale * (B - A)
        backward[i] = B + scale * (A - B)

    if verbose:
        print("Done")

    return forward, backward


if __name__ == "__main__":
    import argparse

    if __name__ == '__main__':

        desc = """
         m    m                      #      m    mm     m mmmmmm
         ##  ##  mmm    m mm  mmmm   # mm    #  #  "m m"      #"
         # ## # #" "#   #"  " #" "#  #"  #    ##    "#"     m#  
         # "" # #   #   #     #   #  #   #   m""m    #     m"   
         #    # "#m#"   #     ##m#"  #   #  m"  "m   #    ##mmmm
                              #                                 
                              " 
         - Interpolate between two XYZ files
         """

        epilog = """
        -----------------------------------------------------
        MINIMAL EXAMPLE
        $ python morph.xyz initial.xyz final.xyz
        
        This will generate an XYZ file containing 
        linearly interpolated geometries, morphed
        from initial to final to initial, both in
        steps of 50.
        -----------------------------------------------------
        DETAILS
        The script first performs an alignment,
        based on the Kabsch algorithm
        (https://en.wikipedia.org/wiki/Kabsch_algorithm).
        The alignment can be slightly tuned by
        the choice of the origin of the optimal
        rotation matrix. By default, center of
        mass is used, but you can choose to use
        the centroids instead.
        
        The verbose flag prints out some key
        objects, such as the covariance and
        optimal rotation matrices from the 
        alignment.
        -----------------------------------------------------
        CAVEATS
        These are currently some limitations or unwanted
        behavior:
        - Start and stop geometries MUST have the same
          number of atoms.
        - Atoms in start and stop should match in a 
          pair-wise manner (line 1 in start should corre-
          spond to line 1 in stop, and so on). If this is
          not the case, then weird atomic transitions may
          take place (e.g., H atom in start becomes C atom
          in stop, but the visualizing software is not able
          to tell, so the representation stays like that 
          for an H atom).
        -----------------------------------------------------
        AUTHOR
        Anders Brakestad
        PhD Candidate in Computational Chemistry
        University of Tromsø The Arctic University of Norway
        """

        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         prog="morph.py",
                                         description=desc,
                                         epilog=epilog)

        parser.add_argument("start", type=str, metavar="<start.xyz>", help="Start geometry of interpolation")
        parser.add_argument("stop", type=str, metavar="<stop.xyz>", help="Final geometry of interpolation")
        parser.add_argument("-f", "--fname", type=str, metavar="str", default="morphed.xyz",
                            help="Filename of output morph xyz file (default: morphed.xyz)")
        parser.add_argument("-d", "--dest", type=str, metavar="str", default=".", help="Output directory (default: .)")
        parser.add_argument("-n", "--nsteps", type=int, default=50, metavar="int",
                            help="Number of interpolation steps (default: 50)")
        parser.add_argument("-o", "--origin", type=str, metavar="str",
                            choices=["centroid", "com"], default="com",
                            help="Origin for alignment {centroid/com} (default: com)")
        parser.add_argument("-l", "--loop", action="store_true",
                            help="Don't loop the morph back to start (default: true)")
        parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output (default: false)")
        args = parser.parse_args()

        # Check file types
        assert all([args.start.split(".")[-1] == "xyz",
                    args.stop.split(".")[-1] == "xyz"]), "Geometries must be in XYZ files"

        # Read in geometries to numpy arrays
        p, lp = get_coords(args.start)
        q, lq = get_coords(args.stop)

        # Check number of atoms
        assert p.shape == q.shape, "Start and stop must have the same number of atoms!"

        # Align structures based on Kabsch algorithm
        p = align(p, q, verbose=args.verbose, method=args.origin, la=lp, lb=lq)

        # Generate the interpolated arrays, morphed to target and back again to probe
        forward, backward = interpolate(p, q, nsteps=args.nsteps, verbose=args.verbose)

        # Write geometries to file
        with open(os.path.join(args.dest, args.fname), "w") as f:
            if args.verbose:
                print("Writing forward interpolation...")
            for step in forward:
                f.write(f"{step.shape[0]}\n")
                f.write(f"Forward morphing\n")
                for label, atom in zip(lp, step):
                    f.write(f"{label} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n")
            if args.loop:
                if args.verbose:
                    print("Writing backward interpolation...")
                for step in backward:
                    f.write(f"{step.shape[0]}\n")
                    f.write(f"Backward morph\n")
                    for label, atom in zip(lp, step):
                        f.write(f"{label} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n")