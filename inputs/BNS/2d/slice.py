"""
'prepare' reads the AthenaK binary slice (dens, temperature, s_00 = Ye), maps
every cell to the thermodynamic quantities needed by bns_nurates (nb, T, Yn,
Yp, mu_e, mu_n, mu_p) by trilinear interpolation of the CompOSE table in
(log nb, Yq, log T), exactly as done by AthenaK, and writes:
    - slice_points.txt:   input of tests_cost/tests_opacities_cost.cpp
    - slice_geometry.npz: meshblock geometry, useful to make 2d plots

Usage:
    python slice.py prepare --snapshot bns.prim_xz.NNNNN.bin
                            --eos DD2.h5 --out-dir DIR
"""

import argparse
import os
import numpy as np

# CGS constants used by AthenaK (src/eos/primitive-solver/unit_system.cpp)
C_CGS = 2.99792458e10  # [cm s^-1]
G_CGS = 6.67408e-8  # [cm^3 g^-1 s^-2]
MSUN_CGS = 1.98848e33  # [g]
MEV_CGS = 1.6021766208e-6  # [erg]


def read_athenak_binary(filename):
    """Read an AthenaK binary output (format version 1.1)"""
    with open(filename, "rb") as fp:
        fp.seek(0, 2)
        filesize = fp.tell()
        fp.seek(0, 0)

        code_header = fp.readline().split()
        if code_header[0] != b"Athena" or code_header[-1].split(b"=")[-1] != b"1.1":
            raise TypeError(f"{filename} is not an AthenaK binary file (v1.1)")

        pheader_count = int(fp.readline().split(b"=")[-1])
        pheader = {}
        for _ in range(pheader_count - 1):
            key, val = [x.strip() for x in fp.readline().decode("utf-8").split("=")]
            pheader[key] = val
        locsize = int(pheader["size of location"])
        varsize = int(pheader["size of variable"])
        locdtype = np.float64 if locsize == 8 else np.float32
        vardtype = np.float64 if varsize == 8 else np.float32

        nvars = int(fp.readline().split(b"=")[-1])
        var_names = [v.decode("utf-8") for v in fp.readline().split()[1:]]
        header_size = int(fp.readline().split(b"=")[-1])
        header = [
            line.decode("utf-8").split("#")[0].strip()
            for line in fp.read(header_size).split(b"\n")
        ]
        nghost = None
        block = None
        for line in header:
            if line.startswith("<"):
                block = line
            elif block == "<mesh>" and line.split("=")[0].strip() == "nghost":
                nghost = int(line.split("=")[1])

        blocks = []
        while fp.tell() < filesize:
            idx = np.frombuffer(fp.read(24), dtype=np.int32).astype(np.int64) - nghost
            fp.read(16)  # logical location and level
            geom = np.frombuffer(fp.read(6 * locsize), dtype=locdtype).astype(np.float64)
            nx1 = idx[1] - idx[0] + 1
            nx2 = idx[3] - idx[2] + 1
            nx3 = idx[5] - idx[4] + 1
            data = np.fromfile(fp, dtype=vardtype, count=nvars * nx1 * nx2 * nx3)
            data = data.reshape(nvars, nx3, nx2, nx1).astype(np.float64)
            blocks.append(
                {
                    "geom": geom,
                    "data": {v: data[i] for i, v in enumerate(var_names)},
                }
            )

    return float(pheader["time"]), blocks


def prepare(args):
    import h5py
    from scipy.interpolate import RegularGridInterpolator

    time, blocks = read_athenak_binary(args.snapshot)
    print(f"Snapshot: {args.snapshot}")
    print(f"Time = {time:.2f} M_sun, {len(blocks)} meshblocks")

    # Slice plane: the two directions with more than one cell
    shape = blocks[0]["data"]["dens"].shape  # (nx3, nx2, nx1)
    if shape[1] == 1:
        plane, ia, ib = "xz", 0, 4  # x1 along columns, x3 along rows
    elif shape[0] == 1:
        plane, ia, ib = "xy", 0, 2
    else:
        raise ValueError("Expected a 2D slice")

    rho, temp, ye = [], [], []
    geometry = []
    offset = 0
    for b in blocks:
        d = b["data"]
        rho.append(np.squeeze(d["dens"]).ravel())
        temp.append(np.squeeze(d["temperature"]).ravel())
        ye.append(np.squeeze(d["s_00"]).ravel())
        nrow, ncol = np.squeeze(d["dens"]).shape
        g = b["geom"]
        geometry.append([g[ia], g[ia + 1], g[ib], g[ib + 1], ncol, nrow, offset])
        offset += nrow * ncol
    rho = np.concatenate(rho)
    temp = np.concatenate(temp)
    ye = np.concatenate(ye)

    # EOS table (CompOSE format, as used by AthenaK)
    with h5py.File(args.eos, "r") as f:
        tab_nb = f["nb"][:]  # [fm^-3]
        tab_yq = f["yq"][:]
        tab_t = f["t"][:]  # [MeV]
        mn = float(f["mn"][()])  # [MeV]
        fields = {k: f[k][:] for k in ["Y[n]", "Y[p]", "Q3", "Q4", "Q5"]}

    # Baryon number density: AthenaK uses mb = mn
    length_unit = G_CGS * MSUN_CGS / C_CGS**2  # [cm]
    rho_cgs = rho * MSUN_CGS / length_unit**3  # [g cm^-3]
    nb_fm = rho_cgs / (mn * MEV_CGS / C_CGS**2) * 1e-39  # [fm^-3]

    # Trilinear interpolation in (log nb, yq, log T), clamping to the table range
    grid = (np.log(tab_nb), tab_yq, np.log(tab_t))
    pts = np.column_stack(
        [
            np.log(np.clip(nb_fm, tab_nb[0], tab_nb[-1])),
            np.clip(ye, tab_yq[0], tab_yq[-1]),
            np.log(np.clip(temp, tab_t[0], tab_t[-1])),
        ]
    )
    interp = {
        k: RegularGridInterpolator(grid, v, method="linear")(pts)
        for k, v in fields.items()
    }

    mu_b = mn * (interp["Q3"] + 1.0)
    mu_q = mn * interp["Q4"]
    mu_l = mn * interp["Q5"]
    mu_n = mu_b
    mu_p = mu_b + mu_q
    mu_e = mu_l - mu_q

    os.makedirs(args.out_dir, exist_ok=True)
    pts_file = os.path.join(args.out_dir, "slice_points.txt")
    out = np.column_stack(
        [
            nb_fm * 1e18,  # [nm^-3]
            temp,
            interp["Y[n]"],
            interp["Y[p]"],
            mu_e,
            mu_n,
            mu_p,
        ]
    )
    header = (
        f"snapshot {os.path.abspath(args.snapshot)}, time = {time:.6e} M_sun\n"
        f"eos {os.path.abspath(args.eos)}\n"
        "nb [nm^-3], T [MeV], Yn, Yp, mu_e [MeV], mu_n [MeV], mu_p [MeV]"
    )
    np.savetxt(pts_file, out, fmt="%.10e", header=header)

    np.savez(
        os.path.join(args.out_dir, "slice_geometry.npz"),
        time=time,
        plane=plane,
        geometry=np.array(geometry),
        rho_cgs=rho_cgs,
        temp=temp,
        ye=ye,
    )

    print(f"Number of points = {len(nb_fm)}")
    print(f"rho range = [{rho_cgs.min():.3e}, {rho_cgs.max():.3e}] g cm^-3")
    print(f"T range   = [{temp.min():.3e}, {temp.max():.3e}] MeV")
    print(f"Ye range  = [{ye.min():.3e}, {ye.max():.3e}]")
    print(f"Written {pts_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_prepare = subparsers.add_parser("prepare", help="Prepare the slice points")
    p_prepare.add_argument("--snapshot", required=True, help="AthenaK binary slice")
    p_prepare.add_argument("--eos", required=True, help="CompOSE EOS table (HDF5)")
    p_prepare.add_argument("--out-dir", default=".", help="Output directory")
    p_prepare.set_defaults(func=prepare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

