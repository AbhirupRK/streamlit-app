#%% This is a Streamlit web application for visualizing atomic magnetic moments from DFT output files.

import io
from pathlib import Path
from modules.dft import ReadVaspOutput, ReadSiestaOutput
import numpy as np
from ase.io import write
from ase import Atoms
import py3Dmol
import streamlit as st


def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()


def get_dft_output(uploaded_file, program):
    if uploaded_file is not None:
        if program is not None:
            try:
                if program.lower() == "vasp":
                    dft_output = ReadVaspOutput(uploaded_file)
                elif program.lower() == "siesta":
                    dft_output = ReadSiestaOutput(uploaded_file)
            except:
                st.error(f"Couldn't find required information in the uploaded file. Please ensure that the {program} output file is correct or report issue to the developer.")
                dft_output = None

        else:
            st.error("Please select the DFT code.")
            dft_output = None
    else:
        st.error("Please upload a DFT output file to start with.")
        st.info("Currently, it only supports VASP and Siesta output files (e.g. OUTCAR, siesta.out). Note that the DFT calculation should be either in spin-polarized or in non-collinear mode.")
        dft_output = None
    return dft_output


def process_atoms(dft_output, species_to_plot):

    species = dft_output.data['species']
    nions = dft_output.data["nions"]
    lattice = dft_output.data["lattice"]
    pos_cart = dft_output.data["pos_cart"]
    magmom = dft_output.data["magmom"]

    structure = Atoms(
        symbols=[sym for sym, count in zip(species, nions) for _ in range(count)],
        positions=pos_cart,
        cell=lattice,
        pbc=True
    )
    if len(species_to_plot)>0:
        atom_indices_to_plot = [i for i, sym in enumerate(structure.symbols) if sym in species_to_plot]
    else:
        st.error("ERROR: Please select at least one species to plot.")
        st.stop()

    return structure[atom_indices_to_plot], magmom[atom_indices_to_plot]


def draw_cell(dft_output):
    # Define the 3 lattice vectors
    origin = np.array([0, 0, 0])
    v1, v2, v3 = dft_output.data["lattice"]

    # Generate the 8 vertices of the parallelepiped
    cell_vertices = np.array([
        origin,
        origin + v1,
        origin + v2,
        origin + v3,
        origin + v1 + v2,
        origin + v1 + v3,
        origin + v2 + v3,
        origin + v1 + v2 + v3
    ])

    # Edges of the cube
    cell_edges = [
        (0, 1), (0, 2), (0, 3),  # Edges from origin
        (1, 4), (1, 5),          # Edges from v1
        (2, 4), (2, 6),          # Edges from v2
        (3, 5), (3, 6),          # Edges from v3
        (4, 7), (5, 7), (6, 7)   # Top edges
    ]
    return cell_vertices, cell_edges


def plot_structure(structure, magmom, cell_vertices, cell_edges):

    buffer = io.StringIO()
    write(buffer, structure, format='xyz')

    view = py3Dmol.view(width='100%', height=700, style={'backgroundColor': 'white'})
    view.addModel(buffer.getvalue(), "xyz", {"doAssembly": True}); buffer.close()
    if plot_bonds:
        view.setStyle({'sphere':{'radius':atom_size}, 'stick':{'radius':atom_size/4}})
    else:
        view.setStyle({'sphere':{'radius':atom_size}})

    if plot_cell:
        colors = ['red', 'green', 'blue'] + ['black'] * (len(cell_edges)-3)  # Color for the unit cell lines
        for i in range(len(cell_edges)):
            start_idx, end_idx = cell_edges[i]
            p1 = cell_vertices[start_idx]
            p2 = cell_vertices[end_idx]
            view.addCylinder({
                'start': {'x': float(p1[0]), 'y': float(p1[1]), 'z': float(p1[2])},
                'end': {'x': float(p2[0]), 'y': float(p2[1]), 'z': float(p2[2])},
                'radius': cell_thickness,   # Thickness of the unit cell lines
                'color': colors[i],         # Color of the edges
                'dashed': False             # Solid lines
            })

    if plot_magmom:
        for i, pos in enumerate(structure.positions):
            if min_magmom < np.linalg.norm(magmom[i]) < max_magmom: 
                arrow_start_coords = {'x': float(pos[0]), 'y': float(pos[1]), 'z': float(pos[2])}
                arrow_end_coords = {
                    'x': float(pos[0] + magmom[i,0] * arrow_scale_factor),
                    'y': float(pos[1] + magmom[i,1] * arrow_scale_factor),
                    'z': float(pos[2] + magmom[i,2] * arrow_scale_factor)
                }
                view.addArrow({
                    'start': arrow_start_coords,
                    'end': arrow_end_coords,
                    'radius': arrow_radius,
                    'color': arrow_color,
                    'alpha': 1 # Transparency
                })

    view.zoomTo()
    view.zoom(zoom_factor)
    view.rotate(rot_x, 'x')
    view.rotate(rot_y, 'y')
    view.rotate(rot_z, 'z')
    view.render()

    t = view.js()
    html_content = t.startjs + t.endjs
    return html_content


##### --- Web App Interface --- #####

st.set_page_config(layout="wide") # Optional: Use wide layout for better viewing

# Defining Left and Right columns for the main layout
main_col1, main_col2 = st.columns([0.6, 0.4]) # Main columns for the plot and inputs
with main_col2:
    st.subheader("File Upload & Controls")
    
    # File uploader and DFT code selection
    file_col1, file_col2 = st.columns(2)
    with file_col1:
        uploaded_file = st.file_uploader("Upload DFT output file", type=None, key="uploaded_file")
    with file_col2:
        program = st.selectbox("Select DFT code", options=["VASP", "Siesta"], index=None, key="program")

    dft_output = get_dft_output(uploaded_file, program)
    if dft_output is not None:
        
        default_values = {
            "species_to_plot": dft_output.data['species'],
            "plot_cell": True,
            "plot_bonds": False,
            "plot_magmom": True,
            "min_magmom": 0.001,
            "max_magmom": 1.0,
            "zoom_factor": 1.0,
            "atom_size": 0.7,
            "cell_thickness": 0.06,
            "rot_x": -85,
            "rot_y": 0,
            "rot_z": 15,
            "arrow_scale_factor": 10,
            "arrow_radius": 0.3,
            "arrow_color": "Blue"
        }

        st.session_state.setdefault("species_to_plot", default_values["species_to_plot"])
        st.session_state.setdefault("plot_cell", default_values["plot_cell"])
        st.session_state.setdefault("plot_bonds", default_values["plot_bonds"])
        st.session_state.setdefault("plot_magmom", default_values["plot_magmom"])
        st.session_state.setdefault("min_magmom", default_values["min_magmom"])
        st.session_state.setdefault("max_magmom", default_values["max_magmom"])
        st.session_state.setdefault("zoom_factor", default_values["zoom_factor"])
        st.session_state.setdefault("atom_size", default_values["atom_size"])
        st.session_state.setdefault("cell_thickness", default_values["cell_thickness"])
        st.session_state.setdefault("rot_x", default_values["rot_x"])
        st.session_state.setdefault("rot_y", default_values["rot_y"])
        st.session_state.setdefault("rot_z", default_values["rot_z"])
        st.session_state.setdefault("arrow_scale_factor", default_values["arrow_scale_factor"])
        st.session_state.setdefault("arrow_radius", default_values["arrow_radius"])
        st.session_state.setdefault("arrow_color", default_values["arrow_color"])
        
        # Reset button and species selection
        with file_col2:
            if st.button("Reset to default", use_container_width=True):
                for k, v in default_values.items():
                    if k != "atom_size":
                        st.session_state[k] = v
                if st.session_state["atom_size"] <= default_values["atom_size"]:
                    st.session_state["atom_size"] = st.session_state["atom_size"] + 0.0001
                else:
                    st.session_state["atom_size"] = st.session_state["atom_size"] - 0.0001
            
            species_to_plot = st.multiselect("Show species", options=dft_output.data['species'], key="species_to_plot")
        
        # Configuration options
        config_col1, config_col2, config_col3 = st.columns([7,7,10])
        with config_col1:
            plot_cell = st.checkbox("Unit cell", key="plot_cell")
        with config_col2:
            plot_bonds = st.checkbox("Bonds", key="plot_bonds")
        with config_col3:
            plot_magmom = st.checkbox("Magnetic moments", key="plot_magmom")
        
        # Magnetic moment range
        magmom_col1, magmom_col2 = st.columns(2)
        with magmom_col1:
            min_magmom = st.number_input("Minimum magnetic moment", min_value=0.0, max_value=10.0, step=0.1, format='%0.4f', key="min_magmom")
        with magmom_col2:
            max_magmom = st.number_input("Maximum magnetic moment", min_value=0.0, max_value=10.0, step=0.1, format='%0.4f', key="max_magmom")

        # Zoom factor and atom size
        atom_col1, atom_col2, atom_col3 = st.columns(3)
        with atom_col1:
            zoom_factor = st.number_input("Zoom factor", min_value=0.0, max_value=10.0, step=0.1, key="zoom_factor")
        with atom_col2:
            atom_size = st.number_input("Atom size", min_value=0.1, max_value=2.0, step=0.1, key="atom_size")
        with atom_col3:
            cell_thickness = st.number_input("Cell thickness", min_value=0.0, max_value=0.5, step=0.01, key="cell_thickness")

        # View angles
        rot_col1, rot_col2, rot_col3 = st.columns(3)
        with rot_col1:
            rot_x = st.number_input("Rotation-X", min_value=-180, max_value=180, step=5, key="rot_x")
        with rot_col2:
            rot_y = st.number_input("Rotation-Y", min_value=-180, max_value=180, step=5, key="rot_y")
        with rot_col3:
            rot_z = st.number_input("Rotation-Z", min_value=-180, max_value=180, step=5, key="rot_z")

        # Arrows
        arr_col1, arr_col2, arr_col3 = st.columns(3)
        with arr_col1:
            arrow_scale_factor = st.number_input("Arrow length", min_value=1, max_value=50, step=5, key="arrow_scale_factor")
        with arr_col2:
            arrow_radius = st.number_input("Arrow radius", min_value=0.05, max_value=0.5, step=0.05, key="arrow_radius")
        with arr_col3:
            arrow_color = st.selectbox("Arrow color", options=["Red", "Green", "Blue", "Black"], key="arrow_color").lower()
        
        # Tips to drag figure position
        st.markdown(":violet-badge[:material/emoji_objects: Use **Ctrl + Click** to drag figure position]")

# Main column for displaying the structure or app description
with main_col1:
    if dft_output is not None:
        final_structure, final_magmom = process_atoms(dft_output, species_to_plot)
        cell_vertices, cell_edges = draw_cell(dft_output)
        py3dmol_html = plot_structure(final_structure, final_magmom, cell_vertices, cell_edges)
        st.components.v1.html(py3dmol_html, height=700, width=None, scrolling=False)
    else:
        st.markdown(read_markdown_file("about.md"), unsafe_allow_html=True)

#%%
