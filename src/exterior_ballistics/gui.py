import tkinter
import ttkbootstrap as ttk
from ttkbootstrap.dialogs import Messagebox
from math import pi


from .rangetable import RangeTable, deg_to_rad, dd_to_dms, rad_to_deg
from . import __version__

root = ttk.Window(themename="cyborg")
root.columnconfigure((0, 1, 2), weight=1)
root.resizable(height=True, width=False)
root.title("Artillery Range " + __version__)
style = ttk.Style()


def check_float(val: float) -> bool:
    try:
        float(val)
        cvt = True
    except ValueError:
        cvt = False
    return cvt or val == ""


reg = root.register(check_float)


def populate(frame, inputs) -> list[tkinter.DoubleVar]:
    variables = []
    for i, (input_name, input_unit, input_default) in enumerate(inputs):
        ttk.Label(frame, text=input_name, anchor="e").grid(row=i, column=0, sticky="nsew")
        var = tkinter.DoubleVar(value=input_default)
        variables.append(var)
        ttk.Entry(frame, width=10, textvariable=var, validate="key", validatecommand=(reg, "%P")).grid(
            row=i, column=1, sticky="nsew"
        )
        ttk.Label(frame, text=input_unit).grid(row=i, column=2, sticky="nsew")
    return variables


frame_configs = {
    "Gun/Shell Parameters": [
        ("Velocity", "m/s", 515.0),
        ("Mass", "kg", 21.76),
        ("Caliber", "mm", 122),
        ("i₄₃", "", 1.0),
    ],
    "Mission Parameters": [
        ("Mount Height ASL", "m", 0.0),
        ("Target Height ASL", "m", 0.0),
        ("Dry Temperature", "°C", 15.0),
        ("Ambient Pressure", "kPa", 100.0),
    ],
    "Computation Parameters": [
        ("Time Step", "s", 1.0),
        ("Range Step", "m", 1000),
        ("Min. Elev.", "°", -5),
        ("Max. Elev.", "°", 65),
    ],
}


def setup_frames(frame_configs) -> list[tkinter.DoubleVar]:
    frames_variables = []
    for i, (label, config) in enumerate(frame_configs.items()):
        label_frame = ttk.LabelFrame(root, text=label)
        label_frame.grid(row=0, column=i, sticky="nsew")
        frames_variables.extend(populate(label_frame, config))

    return frames_variables


frame_variables = setup_frames(frame_configs)
columns = {
    "marker": "",
    "range": "Range (m)",
    "elevation": "Elevation (°)",
    "angle": "Angle (°)",
    "velocity": "Velocity (m/s)",
    "time_to_target": "Time (s)",
    "max_height": "Apex (m)",
}


def generate_trees(tree_names):
    trees = []
    for i, tree_name in enumerate(tree_names):
        row = i + 1
        table_frame = ttk.LabelFrame(root, text=tree_name)
        table_frame.grid(row=row, column=0, columnspan=3, sticky="nsew")
        root.rowconfigure(row, weight=1)

        table_frame.columnconfigure(0, weight=1)
        table_frame.rowconfigure(0, weight=1)

        tree = ttk.Treeview(table_frame, columns=list(columns.keys()), show="headings")
        trees.append(tree)
        for key, value in columns.items():
            tree.column(key, anchor="center", stretch=True, width=0)
            tree.heading(key, text=value)
        tree.grid(row=0, column=0, sticky="nsew")

        horizontal_scroll = ttk.Scrollbar(table_frame, orient="horizontal", command=tree.xview)

        vertical_scroll = ttk.Scrollbar(table_frame, orient="vertical", command=tree.yview)
        tree.configure(xscrollcommand=horizontal_scroll.set, yscrollcommand=vertical_scroll.set)

        horizontal_scroll.grid(row=1, column=0, sticky="nsew")
        vertical_scroll.grid(row=0, column=1, sticky="nsew")

    return trees


trees = generate_trees(("Ground Range Table", "Anti-Aircraft Range Table"))


def compute_range_table():
    try:
        v, m, cal_mm, i_43, gun_h, tgt_h, temp_c, pres_kpa, dt_s, dr_m, theta_min, theta_max = (
            v.get() for v in frame_variables
        )
        range_table = RangeTable(
            gun_alt=gun_h,
            tgt_alt=tgt_h,
            ambient_pressure=pres_kpa * 1e3,
            ambient_temperature=temp_c + 273.15,
            min_elev=deg_to_rad(theta_min),
            max_elev=deg_to_rad(theta_max),
            time_increment=dt_s,
            range_increment=dr_m,
            i_43=i_43,
            velocity=v,
            sectional_density=4 * m / (pi * (cal_mm * 1e-3) ** 2),
        )

        for tree, entries in zip(trees, (range_table.ground_entries, range_table.air_entries)):
            for row in tree.get_children():
                tree.delete(row)

            for entry in entries:

                d, m, s = dd_to_dms(rad_to_deg(entry.gun_elevation))
                tree.insert(
                    "",
                    "end",
                    values=(
                        entry.marker,
                        f"{entry.ground_range:.0f}",
                        f"{d: >3.0f}° {abs(m): >2.0f}'",
                        f"{rad_to_deg(entry.impact_angle): >3.0f}°",
                        f"{entry.impact_velocity:.1f}",
                        f"{entry.time_of_flight:.1f}",
                        f"{entry.apex_height:.0f}",
                    ),
                )
    except Exception as e:
        Messagebox.show_error(message=str(e), title="Exception", parent=root, alert=False)
        raise


button = ttk.Button(root, text="Compute Range Table", command=compute_range_table)
button.grid(row=3, column=0, columnspan=3, sticky="nsew")


def grid_configure_recursive(widget, **kwargs):
    stack = list(widget.winfo_children())
    while stack:
        descendent = stack.pop()
        stack.extend(descendent.winfo_children())
        descendent.grid_configure(**kwargs)


grid_configure_recursive(root, padx=2, pady=2)


def main():
    root.mainloop()


if __name__ == "__main__":

    main()
