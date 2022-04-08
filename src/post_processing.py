import numpy as np
from fenics import SubsetIterator
import os
import csv


def calculate_maximum_volume(f, subdomains, subd_id):
    """Minimum of f over subdomains cells marked with subd_id"""
    V = f.function_space()

    dm = V.dofmap()

    subd_dofs = np.unique(
        np.hstack(
            [dm.cell_dofs(c.index()) for c in SubsetIterator(subdomains, subd_id)]
        )
    )

    return np.max(f.vector().get_local()[subd_dofs])


def write_to_csv(derived_quantities_dict, data):
    """
    Exports data to csv according to parameters in derived_quantities_dict

    Arguments:
    - derived_quantities_dict: dict, contains derived quantities parameters
    - data: list, contains the data to be exported
    Returns:
    - True
    """
    if "file" in derived_quantities_dict.keys():
        file_export = ""
        if "folder" in derived_quantities_dict.keys():
            file_export += derived_quantities_dict["folder"] + "/"
            os.makedirs(os.path.dirname(file_export), exist_ok=True)
        if derived_quantities_dict["file"].endswith(".csv"):
            file_export += derived_quantities_dict["file"]
        else:
            file_export += derived_quantities_dict["file"] + ".csv"
        busy = True
        while busy is True:
            try:
                with open(file_export, "w+") as f:
                    busy = False
                    writer = csv.writer(f, lineterminator="\n")
                    for val in data:
                        writer.writerows([val])
            except OSError as err:
                print("OS error: {0}".format(err))
                print(
                    "The file " + file_export + ".txt might currently be busy."
                    "Please close the application then press any key."
                )
                input()

    return True
