import sys
import re
import numpy as np
import matplotlib as mpl
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import cm

mpl.rcParams["hatch.linewidth"] = 0.5
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = [
    r"\usepackage[utf8]{inputenc}",
    r"\usepackage[russian]{babel}",
    r"\usepackage{amsmath}",
]
for d in Path(".").glob("*"):
    if not d.is_dir():
        continue
    title = d.name + " K"
    field = []
    velocity = []
    velocity_err = []
    rates = []
    for f in sorted(d.glob("*.info"), key=lambda x: float(x.name.split(".")[0])):
        field.append(float(f.name.split(".")[0]))
        velocity.append(
            -float(
                re.search(
                    r"Average velocity: { (\S*), \S*, \S* }", f.read_text()
                ).group(1)
            )
        )
        velocity_err.append(
            float(re.search(r"std: { (\S*), \S*, \S* }", f.read_text()).group(1))
        )
        s = re.search(r"Scattering rates: \[(.*)\]", f.read_text()).group(1)
        rates.append([float(x) for x in s.split(",")])

    field = np.array(field) / 100 / 1000
    velocity = np.array(velocity) * 100 / 1e6
    velocity_err = np.array(velocity_err) * 100 / 1e6
    rates = np.array(rates).T / 1e13
    rates[[0, 1]] = rates[[1, 0]]

    plt.cla()
    plt.xlim(field.min(), field.max())
    plt.ylim(0, velocity.max() * 1.2)
    plt.xlabel(r"Электрическое поле, кВ/см")
    plt.ylabel(r"Скорость дрейфа носителей, $10^{6}\ \text{см}/\text{с}$")
    plt.errorbar(field, velocity, velocity_err, fmt="k.", ecolor="r", capsize=2)
    fit = np.polyfit(field[:3], velocity[:3], 1)
    plt.plot(
        field,
        np.polyval(fit, field),
        "k-",
        lw=1,
        label=r"$\mu = %d\ \text{см}^2 / (\text{В}\cdot\text{с})$"
        % int(fit[0] * 1e6 / 1e3),
    )
    plt.legend(loc="lower right")
    plt.savefig(str(d) + "_vah.png")

    plt.cla()
    plt.xlabel(r"Электрическое поле, кВ/см")
    plt.ylabel(r"Вероятность рассеяния, $10^{13}\ \text{с}^{-1}$")
    stacks = plt.stackplot(field, rates)
    hatches = ["xx", "/////"] + [r"\\\\"] * len(rates[2:])
    handles = stacks[:2] + [stacks[-1]]
    labels = ["Примеси", "Акустические фононы", "Оптические фононы"]
    for s, h in zip(stacks, hatches):
        if h == r"\\\\" and s != stacks[-1]:
            s.set_linewidth(0)
        else:
            s.set_linewidth(0.7)
        s.set_facecolor("none")
        s.set_alpha(1)
        s.set_hatch(h)
    plt.xlim(field.min(), field.max())
    plt.ylim(0, 1.3 * rates.sum(axis=0).max())
    plt.legend(handles, labels)
    plt.savefig(str(d) + "_rates.png")
