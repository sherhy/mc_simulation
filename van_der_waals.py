import shelve

from plotter import scatter_van_der_waals, clear_plt, poly_approx
from numpy import log


def plot_van_der_waals(kts):
    print("creating Van Der Waals.png")
    all_v_val = list()
    all_p_val = list()
    skip_list = ['4.80_5', '4.80_10', '3.54_10', '2.51_10', '1.80_10',
                 '1.39_7', ]
    for kt in kts:
        v_values = list()
        p_values = list()
        with shelve.open(f"./db/vdw_{kt}") as vdw:
            keys = list(vdw.keys())
            keys.sort()
            # print(keys)
            for key in keys:
                if key in skip_list:
                    continue
                v_values.append(vdw[key][0])
                all_v_val.append(vdw[key][0])
                p_values.append(vdw[key][1])
                all_p_val.append(vdw[key][1])
        # print(v_values)
        # print(p_values)
        clear_plt()
        poly_approx(v_values, p_values)
        scatter_van_der_waals(v_values, p_values, f"Van der Waals Isotherms t={kt}")

    clear_plt()
    scatter_van_der_waals(all_v_val, all_p_val, "Van der Waals Isotherms")


def maxwell_construction(x, y, name):
    pass
    # scatter_van_der_waals(all_v_val, all_p_val, "Van der Waals Isotherms")


if __name__ == "__main__":
    pass
    # temps = [0.50]
    # plot_van_der_waals(temps)
