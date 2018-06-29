from pathlib import Path
import hashlib
import yaml
import subprocess
import pickle
import numpy as np

from plot import myPlot

from runcontrols_v2 import do_stuff


def make_results(a: float, b: float, c: float, d: float, end_time: int = 40):
    """Make results."""
    with open("models/my_template.txt", "r") as template_handle:
        template = template_handle.read() % (a, b, c, d)
    with open("models/L5PCbiophys3.hoc", "w") as out_handle:
        out_handle.write(template)

    subprocess.run("/home/jakobes/dev/neuron/nrn/bin/nrnivmodl")
    biophys_file = "models/L5PCbiophys3{}.hoc".format("")   # Want to change this file
    template_file = "models/L5PCtemplate.hoc"
    icell = 0
    morphology_file = "morphologies/cell{}.asc".format(icell + 1)
    picklelist = do_stuff(icell, biophys_file, template_file, morphology_file, end_time)
    return picklelist


def get_data(a: float, b: float, c: float, d: float, end_time: int = 1000):
    """
    Take four parameters and see whether this case has already been computed.

    Arguments:
        a, b, c, d: The conductivities to the NaTa_modt hoc file.

    Returns
    """
    m = hashlib.sha1()
    m.update("{}{}{}{}".format(a, b, c, d).encode())
    my_key = m.hexdigest()[:6]      # probably unique

    model_path = Path("results/model_dict.yaml")

    with open(model_path, "r+") as handle:
        model_dict = yaml.load(handle)
        if my_key in model_dict:
            return load(my_key)

        data = make_results(a, b, c, d, end_time)
        save(my_key, data)
        model_dict[my_key] = a, b, c, d
        yaml.dump(model_dict, handle)
    return data


def save(model_key, data):
    with open("results/control_{}.sav".format(model_key), "wb") as outhandle:
        pickle.dump(data, outhandle)


def load(model_key):
    with open("results/control_{}.sav".format(model_key), "rb") as inhandle:
        return pickle.load(inhandle)


def load_default():
    with open("results/control.sav", "rb") as outhandle:
        return pickle.load(outhandle)


if __name__ == "__main__":
    mutant_fraction = 0.1

    default1 = 2.04 
    default2 = 0.0213
    a = default1*(1 - mutant_fraction)
    b = default1*mutant_fraction

    c = default2*(1 - mutant_fraction)
    d = default2*mutant_fraction

    idx = np.arange(400)
    data = get_data(a, b, c, d)
    time = data[14][0]
    soma = data[15][0]
    dend = data[16][0]

    control = load_default()
    dtime = control[14][0]
    dsoma = control[15][0]
    ddend = control[15][0]

    myPlot(
        (time, soma, "r-", "Soma"),
        (time, dend, "b-", "Dend"),
        (dtime, dsoma, "r--", "Control Soma"),
        (dtime, ddend, "b--", "Control Dend"),
        indices=idx
    )
