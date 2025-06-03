import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

g = ct.Solution('gri30.yaml')
mix = {
    'C2H6': 0.015,
    'C2H4': 0.015,
    'CH4': 0.25,
    'H2': 0.56,
    'CO': 0.08,
    'CO2': 0.03,
    'O2': 0.005,
    'N2': 0.045
}

tot = sum(mix.values())
for sp in mix:
    mix[sp] /= tot

air = {'O2': 0.21, 'N2': 0.79}

T0 = 1100.0
P0 = ct.one_atm

phi_range = np.arange(0.5, 1.51, 0.1)
T_ad = []
dh_mjkg = []
co_vals = []
nox_vals = []
co2_vals = []

print("Equivalence Ratio | Adiabatic Flame Temp [K] | Enthalpy Change [MJ/kg] | CO [ppm] | NOx [ppm] | CO2 [ppm]")
print("------------------------------------------------------------------------------------------------------")

for ph in phi_range:
    try:
        g.set_equivalence_ratio(ph, mix, air)
        g.TP = T0, P0
        h_i = g.enthalpy_mass
        g.equilibrate('HP')
        T_ad.append(g.T)

        co_vals.append(g['CO'].X[0] * 1e6)
        nox_vals.append(sum(g[sp].X[0] for sp in ['NO', 'NO2', 'N2O']) * 1e6)
        co2_vals.append(g['CO2'].X[0] * 1e6)

        g.set_equivalence_ratio(ph, mix, air)
        g.TP = T0, P0
        g.equilibrate('TP')
        h_f = g.enthalpy_mass
        dh = h_i - h_f
        dh_mjkg.append(dh / 1e6)

        print(f"{ph:17.2f} | {g.T:25.2f} | {dh / 1e6:23.6f} | {co_vals[-1]:8.2f} | {nox_vals[-1]:9.2f} | {co2_vals[-1]:9.2f}")

    except Exception as e:
        print(f"{ph:17.2f} | Combustion failed: {str(e)}")
        T_ad.append(np.nan)
        dh_mjkg.append(np.nan)
        co_vals.append(np.nan)
        nox_vals.append(np.nan)
        co2_vals.append(np.nan)

g.TP = 298.15, P0
g.X = mix
h_r = g.enthalpy_mass
g.set_equivalence_ratio(1.0, mix, air)
g.equilibrate('HP')
h_p = g.enthalpy_mass
lhv = abs(h_r - h_p) / 1e6
g.X = mix
g.TP = 273.15, P0
rho = g.density
lhv_vol = lhv * rho

print(f"\nLower Heating Value (LHV): {lhv:.3f} MJ/kg")
print(f"Lower Heating Value (LHV): {lhv_vol:.3f} MJ/nm³")

plt.figure(figsize=(12, 10))
plt.subplot(3, 1, 1)
plt.plot(phi_range, T_ad, marker='o', color='tab:red')
plt.xlabel('Equivalence Ratio (ϕ)')
plt.ylabel('Adiabatic Flame Temperature [K]')
plt.title('Adiabatic Flame Temperature vs. Equivalence Ratio')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.plot(phi_range, co_vals, marker='^', label='CO [ppm]', color='tab:orange')
plt.plot(phi_range, nox_vals, marker='x', label='NOx [ppm]', color='tab:green')
plt.plot(phi_range, co2_vals, marker='d', label='CO2 [ppm]', color='tab:purple')
plt.xlabel('Equivalence Ratio (ϕ)')
plt.ylabel('Emissions [ppm]')
plt.title('Emissions vs. Equivalence Ratio')
plt.legend()
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(phi_range, dh_mjkg, marker='s', color='tab:blue')
plt.xlabel('Equivalence Ratio (ϕ)')
plt.ylabel('Enthalpy Difference [MJ/kg]')
plt.title('Enthalpy Difference vs. Equivalence Ratio')
plt.grid(True)

plt.tight_layout()
plt.show()











