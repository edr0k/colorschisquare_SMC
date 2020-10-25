import os

import numpy as np


def split_file(field, clustername):
    filename = field + '/' + clustername
    filters = ['u', 'f378', 'f395', 'f410', 'f430', 'g', 'f515', 'r', 'f660', 'i', 'f861', 'z']

    with open(filename, 'r') as f:
        # Reads the file with each line as an element of a list
        data = f.readlines()

        # put everything together in a single string
        data = "".join(data)

    # Get the header
    header = data.split(filters[0] + ".fits")[0]
    data = filters[0] + ".fits" + data.split(filters[0] + ".fits")[1]

    # Split and save table for each filter, last filter is done latter
    for i in range(len(filters) - 1):
        # Create new file and save table data
        with open("cluster-phot/" + clustername.upper() + "-" + filters[i], 'w') as f:
            f.write(header)
            f.write(data.split(filters[i + 1] + ".fits")[0])

        # Remove already saved table from data
        data = filters[i + 1] + ".fits" + data.split(filters[i + 1] + ".fits")[1]

    # Save table for the last filter
    with open("cluster-phot/" + clustername.upper() + "-" + filters[-1], 'w') as f:
        f.write(header)
        f.write(data)


param_table = "tables/cat.txt"
cluster_list = np.loadtxt(param_table, usecols=(0), dtype=str)
field_list = np.loadtxt(param_table, usecols=(8), dtype=str)
for i in range(len(cluster_list)):
    path = "fits-splus/" + field_list[i] + "/" + cluster_list[i].lower()
    print(path)
    if (os.path.isfile(path) == True):
        split_file("fits-splus/" + field_list[i], cluster_list[i].lower())
        print(cluster_list)
    # else:
    # print("fail")
