# This file is part of aimsutils.
# (C) 2015 Christoph Schober
import os
import re
from collections import OrderedDict
from aimsutils import find_aimsout

def parse_basis(basisfile="basis-indices.out"):
    """
    Parse the aims basis-indices files to get information on basis functions.

    Parameters
    ----------
    basisfile : str, optional
        The filename for the basis-indices.out file.

    Returns
    -------
    basis_info : OrderedDict
        Dictionary with basis_index: [atom, n, l, m].
    """
    with open(basisfile) as f:
        tmp = f.readlines()[2:]

    all_basis = OrderedDict()
    for basis in tmp:
        idx, basistype, atom, n, l, m = basis.split()
        all_basis[int(idx)] = [int(atom), int(n), int(l), int(m)]
    return all_basis


def parse_aimsout(outfile):
    """
    Parse the aims output file for some information.

    Parameters
    ----------
    outfile : str
        The path to the AIMS output file.

    Returns
    -------
    meta : dict
        Dictionary with restartfile, n_spin, n_empty, n_states, n_basis.
    """
    with open(outfile) as f:
        out = " ".join(f.readlines())
    regex = {"restartfile":
             re.compile(r"Writing restart information to file (\S+)"
                        "| Writing cluster restart information to file (\S+)"),
             "periodic":
             re.compile(r"Found k-point grid: \s+ (\d \s+ \d \s+ \d)"),
             "n_spin":
             re.compile(r"Number of spin channels \s* : \s* ([12]?)"),
             "n_empty":
             re.compile(r"Total number of empty states used during \D* (\d+)"
                        "| Number of empty states per atom: \D* (\d*)"),
             "n_states":
             re.compile(r"Number of Kohn-Sham states \D* (\d*)"),
             "n_states_saved":
             re.compile(r"Restart file: Writing \D* (\d*)"),
             "n_basis":
             re.compile(r"Maximum number of basis functions \D* (\d*)"),
             "proc_gammaK":
             re.compile(r"\| K-points in task\s+(\d):\s+1")}

    meta = dict()
    # Maybe the whole loop should be written explicetly to avoid too many
    # try-except clauses...
    for key, value in regex.items():
        if key == "periodic":
            meta["periodic"] = bool(value.findall(out))
        elif key == "n_empty":
            meta["n_empty"] = int(next(x for x in value.findall(out)[0]
                                   if x != ""))
        elif key == "restartfile":
            meta[key] = [x for x in value.findall(out)[0] if x != ""][0]
        else:
            try:
                meta[key] = int(value.findall(out)[0])
            except ValueError:
                meta[key] = value.findall(out)[0]
            except TypeError:
                pass
            except IndexError:
                # regex found nothing, e.g. periodic in cluster v.v.
                pass

    #meta["n_states_saved"] = min(meta["n_states"],
                                 #meta["n_states"] - meta["n_empty"] + 3)

    if meta["periodic"] is True:
        # assumes aims does 3 numbers for restarts...
        stub = meta["restartfile"][:-3]
        meta["restartfile"] = stub + "{0:03d}".format(meta["proc_gammaK"])

    return meta
