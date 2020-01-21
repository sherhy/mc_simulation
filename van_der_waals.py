import shelve

from plotter import scatter_van_der_waals, clear_plt, poly_approx


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
                p_values.append(vdw[key][1])
        all_v_val.append(v_values)
        all_p_val.append(p_values)
        print(v_values)
        print(p_values)
        clear_plt()
        # poly_approx(v_values, p_values, deg=3)
        scatter_van_der_waals(x=v_values, y=p_values, name=f"Van der Waals Isotherms t={kt}")

    clear_plt()
    for vs, ps in zip(all_v_val, all_p_val):
        poly_approx(xs=vs, ys=ps, deg=3)
        scatter_van_der_waals(x=vs, y=ps, name="Van der Waals Isotherms All", ymax=20)
        # scatter_van_der_waals([],[],name="Van der Waals Isotherms", ymax=10)


def maxwell_construction(x, y, name):
    pass
