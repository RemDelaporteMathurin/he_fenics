from fenics import *
from operator import itemgetter


def subdomains_1D(mesh, materials, size):
    """Iterates through the mesh and mark them
    based on their position in the domain

    Arguments:
        mesh {fenics.Mesh()} -- the mesh
        materials {list} -- contains the dictionaries of the materials
        size {float} -- size of the domain

    Returns:
        fenics.MeshFunction() -- that contains the subdomains
        (0 if no domain was found)
        fenics.MeshFunction() -- that contains the surfaces
        (0 if no domain was found)
    """
    volume_markers = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    for cell in cells(mesh):
        for material in materials:
            if len(materials) == 1:
                volume_markers[cell] = material["id"]
            else:
                if (
                    cell.midpoint().x() >= material["borders"][0]
                    and cell.midpoint().x() <= material["borders"][1]
                ):
                    volume_markers[cell] = material["id"]
    surface_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    surface_markers.set_all(0)
    i = 0
    for f in facets(mesh):
        i += 1
        x0 = f.midpoint()
        surface_markers[f] = 0
        if near(x0.x(), 0):
            surface_markers[f] = 1
        if near(x0.x(), size):
            surface_markers[f] = 2
    return volume_markers, surface_markers


def subdomains(mesh, parameters):
    """Returns two fenics.MeshFunction() for volume and surfaces entities

    Arguments:
        mesh {fenics.Mesh()} -- domain mesh
        parameters {dict} -- contains meshing and materials parameters

    Returns:
        fenics.MeshFunction() -- cell markers
        fenics.MeshFunction() -- facet markers
    """

    mesh_parameters = parameters["mesh_parameters"]
    if "cells_file" in mesh_parameters.keys():
        volume_markers, surface_markers = read_subdomains_from_xdmf(
            mesh, mesh_parameters["cells_file"], mesh_parameters["facets_file"]
        )
    elif "meshfunction_cells" in mesh_parameters.keys() and isinstance(
        mesh_parameters["meshfunction_cells"],
        type(MeshFunction("size_t", mesh, mesh.topology().dim())),
    ):
        volume_markers = mesh_parameters["meshfunction_cells"]
        surface_markers = mesh_parameters["meshfunction_facets"]
    else:
        if "vertices" in mesh_parameters.keys():
            size = max(mesh_parameters["vertices"])
        else:
            size = parameters["mesh_parameters"]["size"]
        if len(parameters["materials"]) > 1:
            check_borders(size, parameters["materials"])
        volume_markers, surface_markers = subdomains_1D(
            mesh, parameters["materials"], size
        )
    return volume_markers, surface_markers


def check_borders(size, materials):
    """Checks that the borders given match

    Arguments:
        size {float} -- float, size of the domain
        materials {list} -- contains dicts with materials parameters

    Raises:
        ValueError: if the borders don't begin at zero
        ValueError: if borders don't match
        ValueError: if borders don't end at size

    Returns:
        bool -- True if everything's alright
    """
    check = True
    all_borders = []
    for m in materials:
        all_borders.append(m["borders"])
    all_borders = sorted(all_borders, key=itemgetter(0))
    if all_borders[0][0] is not 0:
        raise ValueError("Borders don't begin at zero")
    for i in range(0, len(all_borders) - 1):
        if all_borders[i][1] != all_borders[i + 1][0]:
            raise ValueError("Borders don't match to each other")
    if all_borders[len(all_borders) - 1][1] != size:
        raise ValueError("Borders don't match with size")
    return True


def mesh_and_refine(mesh_parameters):
    """Mesh and refine iteratively until meeting the refinement
    conditions.

    Arguments:
        mesh_parameters {dict} -- contains initial number of cells, size,
    and refinements (number of cells and position)

    Returns:
        fenics.Mesh() -- the mesh
    """

    print("Meshing ...")
    initial_number_of_cells = mesh_parameters["initial_number_of_cells"]
    size = mesh_parameters["size"]
    mesh = IntervalMesh(initial_number_of_cells, 0, size)
    if "refinements" in mesh_parameters:
        for refinement in mesh_parameters["refinements"]:
            nb_cells_ref = refinement["cells"]
            refinement_point = refinement["x"]
            print("Mesh size before local refinement is " + str(len(mesh.cells())))
            while len(mesh.cells()) < initial_number_of_cells + nb_cells_ref:
                cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
                cell_markers.set_all(False)
                for cell in cells(mesh):
                    if cell.midpoint().x() < refinement_point:
                        cell_markers[cell] = True
                mesh = refine(mesh, cell_markers)
            print("Mesh size after local refinement is " + str(len(mesh.cells())))
            initial_number_of_cells = len(mesh.cells())
    else:
        print("No refinement parameters found")
    return mesh


def read_subdomains_from_xdmf(mesh, volumetric_file, boundary_file):
    """Reads volume and surface entities from XDMF files

    Arguments:
        mesh {fenics.Mesh()} -- the mesh
        volumetric_file {str} -- path of the XDMF file containing cell
            entities
        boundary_file {str} -- path of the XDMF file containing facet
            entities

    Raises:
        ValueError: if the reading of attribute f in volumetric file fails
        ValueError: if the reading of attribute f in boundary file fails

    Returns:
        fenics.MeshFunction() -- cell markers
        fenics.MeshFunction() -- facet markers
    """

    # Read tags for volume elements
    volume_markers = MeshFunction("size_t", mesh, mesh.topology().dim())
    try:
        XDMFFile(volumetric_file).read(volume_markers, "f")
    except:
        raise ValueError('Attribute should be named "f" in ' + volumetric_file)
    # f is the attribute name carreful

    # Read tags for surface elements
    # (can also be used for applying DirichletBC)
    surface_markers = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    try:
        XDMFFile(boundary_file).read(surface_markers, "f")
    except:
        raise ValueError('Attribute should be named "f" in ' + boundary_file)
    surface_markers = MeshFunction("size_t", mesh, surface_markers)

    print("Succesfully load mesh with " + str(len(volume_markers)) + " cells")
    return volume_markers, surface_markers
