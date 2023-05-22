import math
import os, glob

import matplotlib.pyplot as plt
import numpy as np

time = 1
data_arr = {}
multiplier = 1
T = 0.5

fig, ax = plt.subplots(figsize=(15,15))

def on_press(event):
    global time
    global data_arr
    global multiplier

    if event.key == "e":
        time = time + 1
    if event.key == "w":
        time = time - 1
    ax.clear()
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    for key in data_arr:
        print(key)
        fg = data_arr[key]
        n = fg.grid['N']
        print(fg.get_x())
        ax.plot(fg.get_x(), fg.data[int(time*n/multiplier)])
    ax.set_title(f'Time t= {time*T/multiplier}')
    # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 100 != 0]
    fig.canvas.draw()

fig.canvas.mpl_connect('key_press_event', on_press)


class FuncGrid:

    @staticmethod
    def get_grid(path):
        f = open(path)
        x = f.readline().split(',')
        f.close()
        return {
            'L': float(x[0]),
            'M': int(x[1]),
            'T': float(x[2]),
            'N': int(x[3]),
            'h': float(x[0])/float(x[1]),
            't': float(x[2])/float(x[3])
        }

    def get_x(self):
        return [i*self.grid['h'] for i in range(self.grid['M']+1)]

    def get_t(self):
        return [j*self.grid['t'] for j in range(self.grid['N']+1)]

    def __init__(self, path):
        print("Loading ", path)
        self.grid = FuncGrid.get_grid(path)
        self.data = np.loadtxt(path, delimiter=',', skiprows=1)

def step_up_show():
    global data_arr
    global multiplier
    global T

    T = data_arr['analytic'].grid['T']
    multiplier = 0
    for key in data_arr:
        n = data_arr[key].grid['N']
        if multiplier == 0:
            multiplier = n
        multiplier = math.gcd(multiplier, n)


def calc_dif(fg1, arr_fg):
    n = fg1.grid['N']
    for key in arr_fg:
        n = math.gcd(arr_fg[key].grid['N'], n)

    n_fg1 = fg1.grid['N']
    m_fg1 = fg1.grid['M']
    l1 = {}
    l2 = {}
    linf = {}
    times = []
    for key in arr_fg:
        l1[key] = []
        l2[key] = []
        linf[key] = []

    for i in range(n):
        times.append(i*T/n)

        for key in arr_fg:
            fg = arr_fg[key]
            n_cur = fg.grid['N']
            m_cur = fg.grid['M']
            m = math.gcd(m_cur, m_fg1)

            l1[key].append(np.sum([np.abs(fg.data[int(i*n_cur/n)][int(j*m_cur/m)] - fg1.data[int(i*n_fg1/n)][int(j*m_fg1/m)]) for j in range(m)])/m)
            l2[key].append((np.sum([np.abs(fg.data[int(i*n_cur/n)][int(j*m_cur/m)] - fg1.data[int(i*n_fg1/n)][int(j*m_fg1/m)])**2 for j in range(m)])/m)**(1/2))
            linf[key].append(np.max([np.abs(fg.data[int(i*n_cur/n)][int(j*m_cur/m)] - fg1.data[int(i*n_fg1/n)][int(j*m_fg1/m)]) for j in range(m)]))

    return l1, linf, l2, times

def calc_dif_point(fg1, arr_fg):
    n = fg1.grid['N']
    m = fg1.grid['M']
    for key in arr_fg:
        n = math.gcd(arr_fg[key].grid['N'], n)
        m = math.gcd(arr_fg[key].grid['M'], m)

    n_fg1 = fg1.grid['N']
    m_fg1 = fg1.grid['M']
    dif = {}
    times = []
    pos = []
    for key in arr_fg:
        dif[key] = np.zeros(shape=(n, m))

    for i in range(n):
        times.append(i*T/n)
        for j in range(m):
            pos.append(j*1/m)

            for key in arr_fg:
                fg = arr_fg[key]
                n_cur = fg.grid['N']
                m_cur = fg.grid['M']

                dif[key][i][j] = np.abs(fg.data[int(i*n_cur/n)][int(j*m_cur/m)] - fg1.data[int(i*n_fg1/n)][int(j*m_fg1/m)])

    return dif, times, pos

if __name__ == '__main__':
    # print("Enter path to data file:")
    # path = input()

    anal_path = "data/soliton_up/1e3.0_dataAnalytic.csv"

    lw03_path = "data/soliton_up/1e0.3_dataLW.csv"
    kir03_path = "data/soliton_up/1e0.3_dataKIR.csv"

    lw06_path = "data/soliton_up/1e0.6_dataLW.csv"
    kir06_path = "data/soliton_up/1e0.6_dataKIR.csv"

    lw10_path = "data/soliton_up/1e1.0_dataLW.csv"
    kir10_path = "data/soliton_up/1e1.0_dataKIR.csv"

    lw13_path = "data/soliton_up/1e1.3_dataLW.csv"
    kir13_path = "data/soliton_up/1e1.3_dataKIR.csv"

    lw16_path = "data/soliton_up/1e1.6_dataLW.csv"
    kir16_path = "data/soliton_up/1e1.6_dataKIR.csv"

    lw20_path = "data/soliton_up/1e2.0_dataLW.csv"
    kir20_path = "data/soliton_up/1e2.0_dataKIR.csv"

    lw23_path = "data/soliton_up/1e2.3_dataLW.csv"
    kir23_path = "data/soliton_up/1e2.3_dataKIR.csv"

    lw26_path = "data/soliton_up/1e2.6_dataLW.csv"
    kir26_path = "data/soliton_up/1e2.6_dataKIR.csv"

    lw30_path = "data/soliton_up/1e3.0_dataLW.csv"
    kir30_path = "data/soliton_up/1e3.0_dataKIR.csv"

    lw33_path = "data/soliton_up/1e3.3_dataLW.csv"
    kir33_path = "data/soliton_up/1e3.3_dataKIR.csv"

    lw36_path = "data/soliton_up/1e3.6_dataLW.csv"
    kir36_path = "data/soliton_up/1e3.6_dataKIR.csv"

    data_arr = {
        'analytic': FuncGrid(anal_path),
        'lw03': FuncGrid(lw03_path),
        'kir03': FuncGrid(kir03_path),
        'lw06': FuncGrid(lw06_path),
        'kir06': FuncGrid(kir06_path),
        'lw10': FuncGrid(lw10_path),
        'kir10': FuncGrid(kir10_path),
        'lw13': FuncGrid(lw13_path),
        'kir13': FuncGrid(kir13_path),
        'lw16': FuncGrid(lw16_path),
        'kir16': FuncGrid(kir16_path),
        'lw20': FuncGrid(lw20_path),
        'kir20': FuncGrid(kir20_path),
        'lw23': FuncGrid(lw23_path),
        'kir23': FuncGrid(kir23_path),
        'lw26': FuncGrid(lw26_path),
        'kir26': FuncGrid(kir26_path),
        # 'lw30': FuncGrid(lw30_path),
        # 'kir30': FuncGrid(kir30_path),
        # 'lw33': FuncGrid(lw33_path),
        # 'kir33': FuncGrid(kir33_path),
        # 'lw36': FuncGrid(lw36_path),
        # 'kir36': FuncGrid(kir36_path),
    }

    step_up_show()

    l1_lw, linf_lw, l2_lw, times_lw = calc_dif(data_arr['analytic'], {
        'lw10': data_arr['lw10'],
        'lw13': data_arr['lw13'],
        'lw16': data_arr['lw16'],
        'lw20': data_arr['lw20'],
        'lw23': data_arr['lw23'],
        # 'lw26': data_arr['lw26'],
        # 'lw30': data_arr['lw30'],
        # 'lw33': data_arr['lw33'],
        # 'lw36': data_arr['lw36'],
    })
    l1_kir, linf_kir, l2_kir, times_kir = calc_dif(data_arr['analytic'], {
        'kir10': data_arr['kir10'],
        'kir13': data_arr['kir13'],
        'kir16': data_arr['kir16'],
        'kir20': data_arr['kir20'],
        'kir23': data_arr['kir23'],
        # 'kir26': data_arr['kir26'],
        # 'kir30': data_arr['kir30'],
        # 'kir33': data_arr['kir33'],
        # 'kir36': data_arr['kir36'],
    })

    dif_lw, times_lw1, pos_lw = calc_dif_point(data_arr['analytic'], {
        'lw03': data_arr['lw03'],
        # 'lw06': data_arr['lw06'],
        'lw10': data_arr['lw10'],
        'lw13': data_arr['lw13'],
        'lw16': data_arr['lw16'],
        'lw20': data_arr['lw20'],
        'lw23': data_arr['lw23'],
        # 'lw26': data_arr['lw26'],
        # 'lw30': data_arr['lw30'],
        # 'lw33': data_arr['lw33'],
        # 'lw36': data_arr['lw36'],
    })

    dif_kir, times_kir1, pos_kir = calc_dif_point(data_arr['analytic'], {
        'kir03': data_arr['kir03'],
        # 'kir06': data_arr['kir06'],
        'kir10': data_arr['kir10'],
        'kir13': data_arr['kir13'],
        'kir16': data_arr['kir16'],
        'kir20': data_arr['kir20'],
        'kir23': data_arr['kir23'],
        # 'kir26': data_arr['kir26'],
        # 'kir30': data_arr['kir30'],
        # 'kir33': data_arr['kir33'],
        # 'kir36': data_arr['kir36'],
    })

    # plt.plot(times_lw, np.log(l1_lw['lw10']), 'b')
    # plt.plot(times_lw, np.log(l1_lw['lw13']), 'b')
    # plt.plot(times_lw, np.log(l1_lw['lw16']), 'b')
    # plt.plot(times_lw, np.log(l1_lw['lw20']), 'b')
    #
    # plt.plot(times_kir, np.log(l1_kir['kir10']), 'r')
    # plt.plot(times_kir, np.log(l1_kir['kir13']), 'r')
    # plt.plot(times_kir, np.log(l1_kir['kir16']), 'r')
    # plt.plot(times_kir, np.log(l1_kir['kir20']), 'r')

    pos = 1
    time = 1

    kir_dif = [
        dif_kir['kir03'][time][pos],
        # dif_kir['kir06'][time][pos],
        dif_kir['kir10'][time][pos],
        dif_kir['kir13'][time][pos],
        dif_kir['kir16'][time][pos],
        dif_kir['kir20'][time][pos],
        dif_kir['kir23'][time][pos],
        # dif_kir['kir26'][time][pos],
        # dif_kir['kir30'][time][pos],
        # dif_kir['kir33'][time][pos],
        # dif_kir['kir36'][time][pos]
    ]

    lw_dif = [
        dif_lw['lw03'][time][pos],
        # dif_lw['lw06'][time][pos],
        dif_lw['lw10'][time][pos],
        dif_lw['lw13'][time][pos],
        dif_lw['lw16'][time][pos],
        dif_lw['lw20'][time][pos],
        dif_lw['lw23'][time][pos],
        # dif_lw['lw26'][time][pos],
        # dif_lw['lw30'][time][pos],
        # dif_kir['kir33'][time][pos],
        # dif_kir['kir36'][time][pos]
    ]

    time = 0.1
    # i_lw = np.where(times_lw == 0.1)
    # i_kir = np.where(times_kir == 0.1)
    i_lw = i_kir = 10
    kir_l1 = [
        l1_kir['kir10'][i_kir],
        l1_kir['kir13'][i_kir],
        l1_kir['kir16'][i_kir],
        l1_kir['kir20'][i_kir],
        l1_kir['kir23'][i_kir],
        # l1_kir['kir26'][i_kir],
        # l1_kir['kir30'][i_kir],
        # l1_kir['kir33'][i_kir],
        # l1_kir['kir36'][i_kir],
    ]

    kir_l2 = [
        l2_kir['kir10'][i_kir],
        l2_kir['kir13'][i_kir],
        l2_kir['kir16'][i_kir],
        l2_kir['kir20'][i_kir],
        l2_kir['kir23'][i_kir],
        # l2_kir['kir26'][i_kir],
        # l2_kir['kir30'][i_kir],
        # l2_kir['kir33'][i_kir],
        # l2_kir['kir36'][i_kir],
    ]

    kir_linf = [
        linf_kir['kir10'][i_kir],
        linf_kir['kir13'][i_kir],
        linf_kir['kir16'][i_kir],
        linf_kir['kir20'][i_kir],
        linf_kir['kir23'][i_kir],
        # linf_kir['kir26'][i_kir],
        # linf_kir['kir30'][i_kir],
        # linf_kir['kir33'][i_kir],
        # linf_kir['kir36'][i_kir],
    ]
    #
    lw_l1 = [
        l1_lw['lw10'][i_lw],
        l1_lw['lw13'][i_lw],
        l1_lw['lw16'][i_lw],
        l1_lw['lw20'][i_lw],
        l1_lw['lw23'][i_lw],
        # l1_lw['lw26'][i_lw],
        # l1_lw['lw30'][i_lw],
        # l1_lw['lw33'][i_lw],
        # l1_lw['lw36'][i_lw],
    ]

    lw_l2 = [
        l2_lw['lw10'][i_lw],
        l2_lw['lw13'][i_lw],
        l2_lw['lw16'][i_lw],
        l2_lw['lw20'][i_lw],
        l2_lw['lw23'][i_lw],
        # l2_lw['lw26'][i_lw],
        # l2_lw['lw30'][i_lw],
        # l2_lw['lw33'][i_lw],
        # l2_lw['lw36'][i_lw],
    ]

    lw_linf = [
        linf_lw['lw10'][i_lw],
        linf_lw['lw13'][i_lw],
        linf_lw['lw16'][i_lw],
        linf_lw['lw20'][i_lw],
        linf_lw['lw23'][i_lw],
        # linf_lw['lw26'][i_lw],
        # linf_lw['lw30'][i_lw],
        # linf_lw['lw33'][i_lw],
        # linf_lw['lw36'][i_lw],
    ]
    #
    fineness = [2, 10, 20, 50, 100, 200]
    #
    # print(np.log(fineness), np.log(kir_l1), np.log(lw_l1), np.log(kir_l2), np.log(lw_l2), np.log(kir_linf), np.log(lw_linf))
    print(np.log(fineness), np.log(kir_dif), np.log(lw_dif))

    plt.plot(np.log(fineness), np.log(kir_dif), 'r')
    plt.plot(np.log(fineness), np.log(lw_dif), 'b')
    #
    # plt.plot(np.log(fineness), np.log(kir_l2), 'r^')
    # plt.plot(np.log(fineness), np.log(lw_l2), 'b^')
    #
    # plt.plot(np.log(fineness), np.log(kir_linf), 'rs')
    # plt.plot(np.log(fineness), np.log(lw_linf), 'bs')

    plt.show()