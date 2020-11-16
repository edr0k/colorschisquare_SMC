import matplotlib.pyplot as plt
import numpy as np

dicmod = {'Z': 0, 'Age': 1, 'ebv': 2, 'uJava': 3, 'F378': 4,
          'F395': 5, 'F410': 6, 'F430': 7, 'gSDSS': 8, 'F515': 9, 'rSDSS': 10,
          'F660': 11, 'iSDSS': 12, 'F861': 13, 'zSDSS': 14}

array_Z = (['0001', '0002', '0005', '0010', '0015', '0020', '0030', '0050',
            '0080', '0100', '0152', '0200', '0300', '0400', '0500'])


def load_table(ntable):
    v_table = np.loadtxt(ntable)
    return v_table


def errors(sigma, sizearray):
    yerr = sigma * np.random.randn(sizearray)
    return yerr


def chisquare_core(table, v_obs, v_err, array_band, array_band_ref, dt, pesos):
    # load models or observation
    v_model = table

    # Making the colors of the  model
    array_cor_mod = np.zeros((len(v_model[:, 0]), len(array_band)))
    for i in range(len(array_band)):
        array_cor_mod[:, i] = v_model[:, array_band[i]] - v_model[:, array_band_ref[i]]

    # Making the colors of the data
    array_cor_obs = np.zeros((len(array_band)))
    array_cor_err = np.zeros((len(array_band)))
    array_band = np.array(array_band) - 3
    array_band_ref = np.array(array_band_ref) - 3
    for i in range(len(array_band)):
        array_cor_obs[i] = v_obs[array_band[i]] - v_obs[array_band_ref[i]]
        array_cor_err[i] = np.sqrt(np.square(v_err[array_band[i]]) +
                                   np.square(v_err[array_band_ref[i]]))

    # Matrix of chi^2   
    array_cor_err_max = np.max(array_cor_err)
    array_cor_err_min = np.min(array_cor_err)
    array_peso = pesos
    # print pesos
    # print(array_cor_obs)
    if (array_cor_err_max == 0):
        array_chi = np.sum(array_peso * np.square(array_cor_obs -
                                                  array_cor_mod),
                           axis=1)
    else:
        array_chi = np.sum(array_peso * np.square((array_cor_obs -
                                                   array_cor_mod) / (array_cor_err)),
                           axis=1)

    return array_chi.reshape((len(array_chi), 1))


def cores(v_obs, v_err, array_band, array_band_ref):
    array_band = list(map(dicmod.get, array_band))
    array_band_ref = list(map(dicmod.get, array_band_ref))
    # Making the colors of the data
    array_cor_obs = np.zeros((len(array_band)))
    array_cor_err = np.zeros((len(array_band)))
    array_band = np.array(array_band) - 3
    array_band_ref = np.array(array_band_ref) - 3
    for i in range(len(array_band)):
        array_cor_obs[i] = v_obs[array_band[i]] - v_obs[array_band_ref[i]]
        array_cor_err[i] = np.sqrt(np.square(v_err[array_band[i]]) +
                                   np.square(v_err[array_band_ref[i]]))

    return array_cor_obs, array_cor_err


def chisquare(table, Z_list, v_obs, v_err, bands, band_ref, dt, pesos):
    array_band = list(map(dicmod.get, bands))
    array_band_ref = list(map(dicmod.get, band_ref))

    l_z = len(Z_list)
    l_age = len(table[:, 0:3])
    array_chi = np.column_stack((table[:, 0:3], chisquare_core(table, v_obs, v_err,
                                                               array_band, array_band_ref, dt, pesos)))
    # print array_chi
    pos_min_chi = np.argmin(array_chi[:, 3])
    min_z = array_chi[pos_min_chi, 0]
    min_age = array_chi[pos_min_chi, 1]
    min_ebv = array_chi[pos_min_chi, 2]

    # print("Obtained parameters - log(age)={}, metallicity={}, reddening={}, chi2={}".format(np.log10(min_age),
    #     min_z, min_ebv, array_chi[pos_min_chi, 3]))
    # print("where sun metallicity is 0.0152, {} Z_sun".format(
    #      ((float(min_z)/0.0152))))
    value_age = pos_min_chi - (pos_min_chi / l_age) * l_age
    value_z = str(min_z)[2::]
    while (len(value_z) < 4):
        value_z = value_z + '0'

    return min_z, min_age, min_ebv, array_chi, value_age, value_z


def fit(table, clusters, l, bands, band_ref, pesos):
    Z = clusters[l][1]
    age = clusters[l][2]
    ebv = clusters[l][3]
    dt = 0.05
    sigma = 0
    array_obs = np.array(clusters[l, 3:15])
    array_err = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])
    if (len(clusters[1]) > 15):
        array_err = np.array(clusters[l, 15:27])
    else:
        array_err = np.array([0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02])
    Z_aprox_L = np.argmin(abs(table[:, 0] - Z))
    Z_aprox = table[Z_aprox_L, 0]
    age_aprox_L = np.argmin(abs(table[:, 1] - age))
    age_aprox = table[age_aprox_L, 1]
    ebv_aprox_L = np.argmin(abs(table[:, 2] - ebv))
    ebv_aprox = table[ebv_aprox_L, 2]
    # print(Z_aprox, age_aprox, ebv_aprox)
    for i in range(len(table)):
        if (table[i, 0] == Z_aprox):
            if (table[i, 1] == age_aprox):
                if (table[i, 2] == ebv_aprox):
                    age_L = i
                    # print(age_L)

    array_mod = table[age_L, 3:]
    obs = np.array([array_obs, array_err, array_mod])
    chi = chisquare(table, array_Z, obs[0], obs[1], bands, band_ref, dt, pesos)
    return chi[0], chi[1], chi[2], chi[3]


def full_fit(table1, clusters1, l1, bands1, bands_ref1, pesos1,
             bands2, bands_ref2, pesos2,  # ):
             bands3, bands_ref3, pesos3, age_delta=((0, 0))):
    chi_age = fit(table1, clusters1, l1, bands1, bands_ref1, pesos1)

    table_age = np.array(table1[:, 1])
    table_Z = np.array(table1[:, 0])

    # age_delta = 0
    # age_delta = (0.346*np.log10(chi_age[1]) -2.87) #regressão full
    # age_delta = (0.904*np.log10(chi_age[1]) -8.38) #regressão linear entre 9 e 10
    age_r = np.log10(chi_age[1]) + (age_delta[0] * np.log10(chi_age[1]) + age_delta[1])
    age_aprox_L = np.argmin(abs(np.log10(table1[:, 1]) - age_r))
    age_aprox = np.log10(table1[age_aprox_L, 1])

    # intervalo utilizando idades sem regressão linear
    # age0 = chi_age[1]-10**9
    # agef = chi_age[1]+10**9

    # intervalo utilizando idades com regressão linear para correção de tendência
    age0 = 10 ** (age_aprox - 0.05)
    agef = 10 ** (age_aprox + 0.05)

    # print(np.log10(age0),np.log10(chi_age[1]) - 0.05,np.log10(agef),np.log10(chi_age[1]) + 0.05)
    filter_age = (table_age < agef) & (table_age > age0)
    table2 = table1[filter_age, :]
    # print(np.log10(chi_age[1]),age_error,age_aprox,np.log10(table2[0,1]))
    chi_Z = fit(table2, clusters1, l1, bands2, bands_ref2, pesos2)

    Z0 = chi_Z[0] - 0.001
    Zf = chi_Z[0] + 0.001
    filter_Z = (table_Z > Z0) & (table_Z < Zf)
    table3 = table1[filter_Z, :]
    # print(table3)
    chi_ebv = fit(table3, clusters1, l1, bands3, bands_ref3, pesos3)  # definir bandas e pesos

    return chi_ebv, chi_Z, chi_age


def data_sint(age, Z, ebv, sigma):
    table = load_table('mag.dat')
    Z_aprox_L = np.argmin(abs(table[:, 0] - Z))
    Z_aprox = table[Z_aprox_L, 0]

    age_aprox_L = np.argmin(abs(table[:, 1] - age))
    age_aprox = table[age_aprox_L, 1]

    ebv_aprox_L = np.argmin(abs(table[:, 2] - ebv))
    ebv_aprox = table[ebv_aprox_L, 2]

    for i in range(len(table)):
        if (table[i, 0] == Z_aprox):
            if (table[i, 1] == age_aprox):
                if (table[i, 2] == ebv_aprox):
                    age_L = i
                    # print(age_L)

    array_mod = table[age_L, 3:]

    array_err = errors(sigma, 12)

    array_sint = array_mod + array_err

    array_err = np.ones(len(array_mod)) * sigma

    return array_sint, array_err, array_mod


def make_chi2_matrix(data):
    age_array = np.unique(data[:, 1])
    z_array = np.unique(data[:, 0])

    str_age_array = np.array(['%.2f' % n for n in np.log10(age_array)],
                             dtype=str)
    str_z_array = np.array(['%.4f' % n for n in z_array], dtype=str)
    idx_str_age = np.arange(0, len(age_array) + 1, 5, dtype=int)

    X = np.zeros(len(age_array) * len(z_array)).reshape(len(age_array),
                                                        len(z_array))
    for i in range(len(age_array)):
        for j in range(len(z_array)):
            idx = np.where((data[:, 1] == age_array[i]) &
                           (data[:, 0] == z_array[j]))[0]
            X[i, j] = np.log10(data[idx, 3])

    return X, age_array, z_array, str_age_array, str_z_array, idx_str_age


def multiplot(data):
    table = np.zeros([1065, 4])
    for k in range(0, 21):  # loop avermelhamentos
        for j in range(0, 15):  # Loop metalicidades
            for i in range(0, 71):  # loop idades

                table[i + 71 * j][0] = data[k + 21 * i + 1491 * j][0]
                table[i + 71 * j][1] = data[k + 21 * i + 1491 * j][1]
                table[i + 71 * j][2] = data[k + 21 * i + 1491 * j][2]
                table[i + 71 * j][3] = data[k + 21 * i + 1491 * j][3]

        min_chi = min(data[:, 3])
        max_chi = max(data[:, 3])
        a = make_chi2_matrix(table)

        fig = plt.figure(figsize=(15, 25))
        im = plt.imshow(a[0], interpolation='nearest', vmin=np.log10(min_chi), vmax=np.log10(max_chi),
                        cmap='viridis', aspect=0.1)
        plt.title('Reddening = ' + str(0.01 * k) + '')
        plt.xlabel('metallicity', fontsize=13)
        plt.ylabel('log$_{10}$(age)', fontsize=13)
        plt.xticks(np.arange(len(a[2])), a[4], rotation=0., fontsize=10)
        plt.yticks(a[5], a[3][a[5]], rotation=0., fontsize=10)
        plt.colorbar(im, label=r'$log10(\chi^2$)', cmap='viridis', fraction=0.025)
        plt.show()

# example
# chi=fit(load_table('mag.dat'), np.loadtxt('clusters.txt',usecols=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]), 0,
#       (['uJava', 'F378', 'F395', 'F410', 'F430', 'F515','rSDSS','F660', 'iSDSS', 'F861', 'zSDSS']),
#       (['gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS', 'gSDSS']),
#       np.ones((11)))
# print(chi)
