import pickle
import numpy as np

mydtype = [
    ("current", "f8"),
    ("frac75", "f8"),
    ("frac80", "f8"),
    ("frac85", "f8"),
    ("frac90", "f8"),
    ("frac95", "f8"),
    ("frac100", "f8")
]

data_array = None

prev_current = None
for i in range(6):
    frac = 75 + i*5
    with open("Q1481K_fraction{:d}_percent_f-I_curve.sav".format(frac), "rb") as filehandle:
        data = pickle.load(filehandle)
        current = data[19]
        if prev_current is None:
            prev_current = current
        assert np.array_equal(current, prev_current)
        prev_current = current

        if data_array is None:
            data_array = np.empty(shape=len(current), dtype=mydtype)
            data_array["current"] = current
        data_array["frac{:d}".format(frac)] = data[0][0]


np.save("fI_data.npy", data_array)
