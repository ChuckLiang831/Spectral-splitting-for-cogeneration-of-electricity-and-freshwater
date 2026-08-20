import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



# ------------variables & constants---------
q_sun = 1000  # Wm-2
Ta = 25+273.15  # K
M = 0.018  # kg/mol
h_fg = 2394000  # J/kg
Da = 3e-5  # m2/s
R = 8.314  # J/(molK)
sigma = 5.670374419e-8  # Stefan-Boltzmann constant in W/(m^2K^4)
T_sky = 293.15  # K
epsilon = 0.03  #
PV_len = 7e-2  # m
PV_wide = 5e-2  # m
TA_len = 5e-2  # m
TA_wide = 5e-2  # m
t = 0.01

k = 0.05
ka = 0.026
ha_pv= 4
ha = 4
R_pvside = (1 / ha_pv)
R_TAside = (1/ha)+(t/k)
Tref_PV = 25+273.15
beta = -0.07882
E_PV=0.2
T_PV_down = Ta
T_PV_up = 373.15
b=0.005

def TAmodule(n,q_inlet):
    T_b = [i for i in range(n)]
    T_f = [i for i in range(n)]
    q_cond = [i for i in range(n)]
    q_conv = [0 for _ in range(n)]
    q_rad = [0 for _ in range(n)]
    q_eva = [i for i in range(n)]
    q_side = [i for i in range(n)]
    q_out = [i for i in range(n)]
    q_in = [i for i in range(n)]
    q_net = [i for i in range(n)]
    J_eva = [i for i in range(n)]
    q_guess = [i for i in range(n)]
    q_guess[n - 1] = q_inlet
    T_b[n - 1] = 22+273

    cycle = 1

    while True:
        cycle = cycle + 1
        for i in range(n, 0, -1):
            T_high = 1000
            T_low = Ta
            for iter in range(2000):
                T_mid = (T_high + T_low) / 2
                Psatf = 610.78 * math.exp((17.27 * (T_mid - 273.15)) / (T_mid + 237.3 - 273.15)) #  Pa
                c_fi = Psatf / (R * T_mid)
                if i == n:
                    Psatb = 610.78 * math.exp((17.27 * (Ta - 273.15)) / (Ta + 237.3 - 273.15))  # Pa
                    c_bi = Psatb / (R * Ta)  # Psatb in Pa cf,cb in mol/m3
                    q_condi = ka * (T_mid - Ta) / b
                    q_sidei = ((T_mid + Ta) / 2 - T_sky) / R_TAside
                else:
                    Psatb = 610.78 * math.exp((17.27 * (T_b[i - 1] - 273.15)) / (T_b[i - 1] + 237.3 - 273.15))  # Pa
                    c_bi = Psatb / (R * T_b[i - 1])
                    q_condi = ka * (T_mid - T_b[i - 1]) / b
                    q_sidei = ((T_mid + T_b[i - 1]) / 2 - T_sky) / R_TAside
                J_evai = Da * (c_fi - c_bi) / b # J in mol/m2/s
                q_evai = M * h_fg * J_evai
                q_outi = (TA_len * TA_wide * (q_condi + q_evai) - 4 * TA_len * b * q_sidei) / (TA_len * TA_wide)
                Q_outi = q_outi * TA_len * TA_wide
                q_neti = (Q_outi - TA_len*TA_wide*q_condi) / (TA_wide * TA_len)
                if i == 1:
                    q_radi = sigma * epsilon * (T_mid ** 4 - T_sky ** 4)
                    q_convi = ha * (T_mid - T_sky)
                    q_ini = (TA_len * TA_wide * (q_radi + q_convi + q_outi) + 4 * TA_len * b * q_sidei) / (
                            TA_len * TA_wide)
                else:
                    q_ini = (TA_len * TA_wide * (q_condi + q_evai)) / (TA_len * TA_wide)

                if abs(q_outi - q_guess[i - 1]) < 10e-5:
                    T_low = T_mid
                    q_out[i - 1] = q_outi
                    T_f[i - 1] = T_mid
                    q_in[i - 1] = q_ini
                    q_net[i - 1] = q_neti
                    q_cond[i - 1] = q_condi
                    q_eva[i - 1] = q_evai
                    q_side[i - 1] = q_sidei
                    J_eva[i - 1] = J_evai
                    if i == 1:
                        q_rad[0] = q_radi
                        q_conv[0] = q_convi
                    elif i > 1:
                        q_guess[i - 2] = q_ini
                        q_out[i - 2] = q_ini
                        T_b[i - 2] = T_mid
                    break
                elif q_outi - q_guess[i - 1] > 0:
                    T_high = T_mid
                else:
                    T_low = T_mid
        if abs(q_in[0] - q_inlet) < 10e-5:
            break
        elif cycle > 1000:
            break
        gama = q_inlet / q_in[0]
        q_guess[-1] = q_guess[-1] * gama

    T_b = [x - 273.15 for x in T_b]
    T_f = [x - 273.15 for x in T_f]

    return T_f, T_b, q_cond, q_conv, q_rad, q_eva, q_side, q_out, q_in, q_net, 0, 0, J_eva
def PVmodule(n,q_inlet,eta_ref,beta):
    T_b = [i for i in range(n)]
    T_f = [i for i in range(n)]
    q_cond = [i for i in range(n)]
    q_conv = [0 for _ in range(n)]
    q_rad = [0 for _ in range(n)]
    q_eva = [i for i in range(n)]
    q_side = [i for i in range(n)]
    q_out = [i for i in range(n)]
    q_in = [i for i in range(n)]
    q_net = [i for i in range(n)]
    J_eva = [i for i in range(n)]
    q_guess = [i for i in range(n)]

    T_b[n - 1] = Ta
    T_PV_down = 273.15
    T_PV_up = 373
    i_pv = 1
    while (i_pv<2000):

        for i_Tpv in range(2000):
            T_PV = 0.5 * (T_PV_down + T_PV_up)
            p_PV = 1000 * eta_ref*(1-beta * (T_PV-Tref_PV))   # W / m2；
            qrad_PV = E_PV * sigma * (T_PV ** 4 - Ta ** 4)
            qconv_PV = ha_pv * (T_PV - Ta)
            qloss_PV = (qrad_PV + qconv_PV)  #  W / m2
            q_heatPV = q_inlet - p_PV - qloss_PV  # W / m2

            if q_heatPV < 0:
                T_PV_up=T_PV
            else:
                break

        cycle = 1
        q_guess[n - 1] = 0.1 * q_heatPV
        while True:
            for i in range(n, 0, -1):
                T_high = 373
                T_low = Ta
                for iter in range(2000):
                    T_mid = (T_high + T_low) / 2
                    Psatf = 610.78 * math.exp((17.27 * (T_mid - 273.15)) / (T_mid + 237.3 - 273.15))  # Pa
                    c_fi = Psatf / (R * T_mid)
                    if i == n:
                        Psatb = 610.78 * math.exp((17.27 * (Ta - 273.15)) / (Ta + 237.3 - 273.15))  # Pa
                        c_bi = Psatb / (R * Ta)  # Pa
                        q_condi = ka * (T_mid - Ta) / b
                        q_sidei = ((T_mid + Ta) / 2 - T_sky) / R_pvside
                    else:
                        Psatb = 610.78 * math.exp((17.27 * (T_b[i - 1] - 273.15)) / (T_b[i - 1] + 237.3 - 273.15))  # Pa
                        c_bi = Psatb / (R * T_b[i - 1])
                        q_condi = ka * (T_mid - T_b[i - 1]) / b
                        q_sidei = ((T_mid + T_b[i - 1]) / 2 - Ta) / R_pvside
                    J_evai = Da * (c_fi - c_bi) / b
                    q_evai = M * h_fg * J_evai
                    q_outi = (PV_len * PV_wide * (q_condi + q_evai) - 4 * PV_len * b * q_sidei) / (PV_len * PV_wide)
                    Q_outi = q_outi * PV_len * PV_wide
                    q_neti = (Q_outi - PV_len*PV_wide*q_condi) / (PV_wide * PV_len)
                    if i == 1:
                        q_radi = E_PV * sigma * (T_mid ** 4 - T_sky ** 4)
                        q_convi = ha_pv * (T_mid - Ta)
                        q_ini = (PV_len * PV_wide * (q_radi + q_convi + q_outi) + 4 * PV_len * b * q_sidei) / (
                                PV_len * PV_wide)
                    else:
                        q_ini = (PV_len * PV_wide * (q_condi + q_evai)) / (PV_len * PV_wide)

                    if abs(q_outi - q_guess[i - 1]) < 10e-5:
                        T_low = T_mid
                        q_out[i - 1] = q_outi
                        T_f[i - 1] = T_mid
                        q_in[i - 1] = q_ini
                        q_net[i - 1] = q_neti
                        q_cond[i - 1] = q_condi
                        q_eva[i - 1] = q_evai
                        q_side[i - 1] = q_sidei
                        J_eva[i - 1] = J_evai
                        if i == 1:
                            q_rad[0] = q_radi
                            q_conv[0] = q_convi
                        elif i > 1:
                            q_guess[i - 2] = q_ini
                            q_out[i - 2] = q_ini
                            T_b[i - 2] = T_mid
                        break
                    elif q_outi - q_guess[i - 1] > 0:
                        T_high = T_mid
                    else:
                        T_low = T_mid
            if abs(q_in[0] - q_heatPV) < 10e-4:
                break
            elif cycle > 1000:
                break
            gama = q_heatPV / q_in[0]
            q_guess[-1] = q_guess[-1] * gama
            cycle = cycle + 1

        if abs(T_f[0] - T_PV)<10e-3:
            break
        elif T_f[0] - T_PV < 0:
            T_PV_up = T_PV
            i_pv = i_pv + 1
        elif T_f[0] - T_PV > 0:
            T_PV_down = T_PV
            i_pv = i_pv + 1

    T_f = [x-273.15 for x in T_f ]
    T_b = [x-273.15 for x in T_b ]

    return T_f, T_b, q_cond, qconv_PV, qrad_PV, q_eva, q_side, q_out, q_in, q_net, 0, 0, J_eva, p_PV
    
os.makedirs('./results', exist_ok=True)

# ---------------------------- Fig2. b&c ---------------------------------------------
q = 1000
n_values = [1, 5, 10]
Eg_values = np.arange(0.5, 2.1, 0.1)
eta_values = [0.1394, 0.1806, 0.2281, 0.2513, 0.2864, 0.3082, 0.3223, 0.3274, 0.3257, 0.3291, 0.3164, 0.3014,
              0.2864, 0.2686, 0.2475, 0.2262]
beta_values = [0.00353504, 0.002870889, 0.002409379, 0.00208642, 0.001852307, 0.00167172, 0.001523727, 0.001401779,
               0.001313715, 0.001281761, 0.001342526, 0.001547006, 0.001960585, 0.00266303, 0.003748494, 0.005325519]
q_values = [1001.00, 996.26, 994.38, 957.09, 937.17, 897.63, 858.67, 809.74, 757.99, 718.62, 662.10, 607.52, 557.28,
            508.71, 458.32, 412.83]
ECD_data = {'Eg': Eg_values}

for n in n_values:
    PVther_list = []
    PCE_list = []
    PVT_list = []

    for Eg, eta, beta in zip(Eg_values, eta_values, beta_values):
        result = PVmodule(n, q, eta, beta)
        PVT = result[0][0]
        PVther = sum(result[5]) / 10
        PCE = result[13] / 10

        PVther_list.append(PVther)
        PCE_list.append(PCE)
        PVT_list.append(PVT)

    ECD_data[f'{n}S-ECD eta_ther'] = PVther_list
    ECD_data[f'{n}S-ECD PCE'] = PCE_list
    ECD_data[f'{n}S-ECD T_PV'] = PVT_list

ECD_df = pd.DataFrame(ECD_data)
ECD_df.to_excel('./results/ECD_Eg.xlsx', index=False)

Eg_values = np.arange(0.5, 2.1, 0.1)
eta_values = [
    0.1394, 0.1806, 0.2281, 0.2513, 0.2864, 0.3082, 0.3223, 0.3274, 0.3257, 0.3291,
    0.3164, 0.3014, 0.2864, 0.2686, 0.2475, 0.2262
]
beta_values = [
    0.00353504, 0.002870889, 0.002409379, 0.00208642, 0.001852307, 0.00167172,
    0.001523727, 0.001401779, 0.001313715, 0.001281761, 0.001342526, 0.001547006,
    0.001960585, 0.00266303, 0.003748494, 0.005325519
]
q_values = [
    999.90, 996.26, 994.38, 957.09, 937.17, 897.63, 858.67, 809.74, 757.99, 718.62,
    662.10, 607.52, 557.28, 508.71, 458.32, 412.83
]

SSCD_data = {'Eg': Eg_values}

n_values = [1, 5, 10]
for n in n_values:
    SSCD_PVther_list = []
    SSCD_PCE_list = []
    SSCD_PVT_list = []
    SSCD_TAther_list = []
    SSCD_Total_list = []

    for Eg, eta, beta, q in zip(Eg_values, eta_values, beta_values, q_values):
        result_pv = PVmodule(n, q, eta, beta)
        SSCD_PVT = result_pv[0][0]
        SSCD_PVther = sum(result_pv[5]) / 10
        SSCD_PCE = result_pv[13] / 10

        q2 = 1000 - q
        if n == 5:
            ta_n = 4
        else:
            ta_n = n
        result_ta = TAmodule(ta_n, q2)
        SSCD_TAther = sum(result_ta[5]) / 10
        SSCD_Total = SSCD_PVther + SSCD_TAther

        SSCD_PVther_list.append(SSCD_PVther)
        SSCD_PCE_list.append(SSCD_PCE)
        SSCD_PVT_list.append(SSCD_PVT)
        SSCD_TAther_list.append(SSCD_TAther)
        SSCD_Total_list.append(SSCD_Total)

    SSCD_data[f'{n}S-SSCD eta_PVther'] = SSCD_PVther_list
    SSCD_data[f'{n}S-SSCD PCE'] = SSCD_PCE_list
    SSCD_data[f'{n}S-SSCD T_PV'] = SSCD_PVT_list
    SSCD_data[f'{n}S-SSCD eta_TAther'] = SSCD_TAther_list
    SSCD_data[f'{n}S-SSCD eta_Total'] = SSCD_Total_list

SSCD_df = pd.DataFrame(SSCD_data)
SSCD_df.to_excel('./results/SSCD_Eg.xlsx', index=False)

eta_values = [x * 100 for x in [
    0.1394, 0.1806, 0.2281, 0.2513, 0.2864, 0.3082, 0.3223, 0.3274, 0.3257, 0.3291,
    0.3164, 0.3014, 0.2864, 0.2686, 0.2475, 0.2262
]]
plt.figure(figsize=(14, 10))

plt.subplot(1, 2, 1)
plt.plot(Eg_values, eta_values, color='#154B9A', label='S-Q limit')
plt.plot(Eg_values, ECD_data['1S-ECD PCE'], color='#78A2D5', linestyle='--', label='1S-ECD')
plt.plot(Eg_values, ECD_data['5S-ECD PCE'], color='#F8C165', linestyle='--', label='5S-ECD')
plt.plot(Eg_values, ECD_data['10S-ECD PCE'], color='#E7A5B0', linestyle='--', label='10S-ECD')
plt.plot(Eg_values, SSCD_data['1S-SSCD PCE'], color='#78A2D5', label='1,1S-SSCD')
plt.plot(Eg_values, SSCD_data['5S-SSCD PCE'], color='#F8C165', label='5,4S-SSCD')
plt.plot(Eg_values, SSCD_data['1S-SSCD PCE'], color='#E7A5B0', label='10,10S-SSCD')
plt.xlabel('Eg (eV)')
plt.ylabel('PCE (%)')
plt.ylim(5, 45)
plt.title('Fig. 2b')
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
ax1 = plt.gca()
ax2 = ax1.twinx()

ax1.plot(Eg_values, ECD_data['1S-ECD eta_ther'], color='#C1D4E5', linestyle='--', label='1S-ECD')
ax1.plot(Eg_values, ECD_data['5S-ECD eta_ther'], color='#78A2D5', linestyle='--', label='5S-ECD')
ax1.plot(Eg_values, ECD_data['10S-ECD eta_ther'], color='#4E7DA9', linestyle='--', label='10S-ECD')

ax1.plot(Eg_values, SSCD_data['1S-SSCD eta_Total'], color='#C1D4E5', label='1,1S-SSCD')
ax1.plot(Eg_values, SSCD_data['5S-SSCD eta_Total'], color='#78A2D5', label='5,4S-SSCD')
ax1.plot(Eg_values, SSCD_data['10S-SSCD eta_Total'], color='#4E7DA9', label='10,10S-SSCD')

ax2.plot(Eg_values, ECD_data['1S-ECD T_PV'], color='#E7A5B0', linestyle='--', label='1S-ECD')
ax2.plot(Eg_values, ECD_data['5S-ECD T_PV'], color='#D2527F', linestyle='--', label='5S-ECD')
ax2.plot(Eg_values, ECD_data['10S-ECD T_PV'], color='r', linestyle='--', label='10S-ECD')
ax2.plot(Eg_values, SSCD_data['1S-SSCD T_PV'], color='#E7A5B0', label='1,1S-SSCD')
ax2.plot(Eg_values, SSCD_data['5S-SSCD T_PV'], color='#D2527F', label='5,4S-SSCD')
ax2.plot(Eg_values, SSCD_data['10S-SSCD T_PV'], color='r', label='10,10S-SSCD')

ax1.set_xlabel('Eg (eV)')
ax1.set_ylabel(r'$\eta_{ther}$ (%)')
ax2.set_ylabel(r'$T_{PV}$ (°C)')
ax1.set_ylim(40, 300)
ax2.set_ylim(20, 100)
plt.title('Fig. 2c')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.grid(True)
plt.tight_layout()
plt.savefig('./results/Fig. 2 b&c.png')

# ------------------------Fig.2 d&e-------------------------
q = 1000
eta_values = [0.21]
beta_values = [0.00185231]
b_values = np.arange(0.001, 0.009 + 0.0005, 0.0005)
ECD_data2 = {'b': b_values}

n_values = [1, 5, 10]
for n in n_values:
    PVther_list = []
    PCE_list = []
    PV_T_list = []

    for eta, beta in zip(eta_values, beta_values):
        for b in b_values:
            result_pv = PVmodule(n, q, eta, beta)
            PV_T = result_pv[0][0]

            PVther = sum(result_pv[5]) / 10

            PCE = result_pv[13] / 10

            PV_T_list.append(PV_T)
            PVther_list.append(PVther)
            PCE_list.append(PCE)

    ECD_data2[f'{n}S-ECD T_PV'] = PV_T_list
    ECD_data2[f'{n}S-ECD eta_ther'] = PVther_list
    ECD_data2[f'{n}S-ECD PCE'] = PCE_list

ECD_df = pd.DataFrame(ECD_data2)
ECD_df.to_excel('./results/ECD_b.xlsx', index=False)

q = 697
eta_values = [0.21]
beta_values = [0.00185231]
b_values = np.arange(0.001, 0.009 + 0.0005, 0.0005)
SSCD_data2 = {'b': b_values}

n_values = [1, 5, 10]
for n in n_values:
    PVther_list = []
    PCE_list = []
    PV_T_list = []
    Total_list = []

    for eta, beta in zip(eta_values, beta_values):
        for b in b_values:
            result_pv = PVmodule(n, q, eta, beta)
            PV_T = result_pv[0][0]
            PVther = sum(result_pv[5]) / 10
            PCE = result_pv[13] / 10
            q2 = 1000 - q
            if n == 5:
                ta_n = 4
            else:
                ta_n = n

            result_ta = TAmodule(ta_n, q2)
            TAther = sum(result_ta[5]) / 10
            Total = PVther + TAther

            PV_T_list.append(PV_T)
            PVther_list.append(PVther)
            PCE_list.append(PCE)
            Total_list.append(Total)

    SSCD_data2[f'{n}S-SSCD T_PV'] = PV_T_list
    SSCD_data2[f'{n}S-SSCD eta PVther'] = PVther_list
    SSCD_data2[f'{n}S-SSCD PCE'] = PCE_list
    SSCD_data2[f'{n}S-SSCD eta_Total'] = Total_list

SSCD_df = pd.DataFrame(SSCD_data2)
SSCD_df.to_excel('./results/SSCD_b.xlsx', index=False)

plt.figure(figsize=(14, 10))

plt.subplot(1, 2, 1)
ax1 = plt.gca()
ax2 = ax1.twinx()

ax1.plot(b_values, ECD_data2['1S-ECD PCE'], color='#78A2D5', marker='s', markerfacecolor='none',
         markeredgecolor='#78A2D5', linestyle='-', label='1S-ECD')
ax1.plot(b_values, ECD_data2['5S-ECD PCE'], color='#F8C165', marker='s', markerfacecolor='none',
         markeredgecolor='#F8C165', linestyle='-', label='5S-ECD')
ax1.plot(b_values, ECD_data2['10S-ECD PCE'], color='#E7A5B0', marker='s', markerfacecolor='none',
         markeredgecolor='#E7A5B0', linestyle='-', label='10S-ECD')
ax1.plot(b_values, SSCD_data2['1S-SSCD PCE'], color='#78A2D5', marker='o', markerfacecolor='none',
         markeredgecolor='#78A2D5', linestyle='-', label='1,1S-SSCD')
ax1.plot(b_values, SSCD_data2['5S-SSCD PCE'], color='#F8C165', marker='o', markerfacecolor='none',
         markeredgecolor='#F8C165', linestyle='-', label='5,4S-SSCD')
ax1.plot(b_values, SSCD_data2['10S-SSCD PCE'], color='#E7A5B0', marker='o', markerfacecolor='none',
         markeredgecolor='#E7A5B0', linestyle='-', label='10,10S-SSCD')

ax2.plot(b_values, ECD_data2['1S-ECD T_PV'], color='#78A2D5', marker='s', markerfacecolor='#78A2D5',
         linestyle='-', label='1S-ECD')
ax2.plot(b_values, ECD_data2['5S-ECD T_PV'], color='#F8C165', marker='s', markerfacecolor='#F8C165',
         linestyle='-', label='5S-ECD')
ax2.plot(b_values, ECD_data2['10S-ECD T_PV'], color='#E7A5B0', marker='s', markerfacecolor='#E7A5B0',
         linestyle='-', label='10S-ECD')
ax2.plot(b_values, SSCD_data2['1S-SSCD T_PV'], color='#78A2D5', marker='o', markerfacecolor='#78A2D5',
         linestyle='-', label='1,1S-SSCD')
ax2.plot(b_values, SSCD_data2['5S-SSCD T_PV'], color='#F8C165', marker='o', markerfacecolor='#F8C165',
         linestyle='-', label='5,4S-SSCD')
ax2.plot(b_values, SSCD_data2['10S-SSCD T_PV'], color='#E7A5B0', marker='o', markerfacecolor='#E7A5B0',
         linestyle='-', label='10,10S-SSCD')

ax1.set_xlabel('b (m)')
ax1.set_ylabel('PCE (%)')
ax1.set_ylim(19, 25)
ax2.set_ylabel(r'$T_{PV}$ (°C)')
ax2.set_ylim(0, 80)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.subplot(1, 2, 2)
plt.plot(b_values, ECD_data2['1S-ECD eta_ther'], color='#78A2D5', marker='s', markerfacecolor='none',
         markeredgecolor='#78A2D5', linestyle='-', label='1S-ECD')
plt.plot(b_values, ECD_data2['5S-ECD eta_ther'], color='#F8C165', marker='s', markerfacecolor='none',
         markeredgecolor='#F8C165', linestyle='-', label='5S-ECD')
plt.plot(b_values, ECD_data2['10S-ECD eta_ther'], color='#E7A5B0', marker='s', markerfacecolor='none',
         markeredgecolor='#E7A5B0', linestyle='-', label='10S-ECD')
plt.plot(b_values, SSCD_data2['1S-SSCD eta_Total'], color='#78A2D5', marker='o', markerfacecolor='none',
         markeredgecolor='#78A2D5', linestyle='-', label='1,1S-SSCD')
plt.plot(b_values, SSCD_data2['5S-SSCD eta_Total'], color='#F8C165', marker='o', markerfacecolor='none',
         markeredgecolor='#F8C165', linestyle='-', label='5,4S-SSCD')
plt.plot(b_values, SSCD_data2['10S-SSCD eta_Total'], color='#E7A5B0', marker='o', markerfacecolor='none',
         markeredgecolor='#E7A5B0', linestyle='-', label='10,10S-SSCD')

plt.xlabel('b (m)')
plt.ylabel(r'$\eta_{ther}$ (%)')
plt.ylim(0, 500)
plt.legend()
plt.tight_layout()
plt.savefig('./results/Fig. 2 d&e.png')
