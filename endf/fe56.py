import numpy as np
import matplotlib.pyplot as plt
from ENDFtk import sequence

# Caminho do arquivo ENDF
path = "/mnt/medley/LucasAnalysis/neutrons-version.VIII.1/n-026_Fe_056.endf"

# Carrega o arquivo diretamente
tape = sequence.read(path)

mat_number = 26056
material = None

for mat in tape:
    if mat.MAT == mat_number:
        material = mat
        break

if material is None:
    raise RuntimeError(f"Material MAT={mat_number} não encontrado no arquivo.")

mf3 = material.MF(3)

if not mf3.hasMT(102):
    raise RuntimeError("MT=102 (n,p) não encontrado no MF=3 do material.")

mt102 = mf3.MT(102)

data = mt102.section(1).data()

energies = np.array(data.x())
xs = np.array(data.y())

mask = (energies >= 1.0) & (energies <= 40.0)
energies = energies[mask]
xs = xs[mask]

plt.figure(figsize=(8,6))
plt.plot(energies, xs, label='Fe-56 (n,xp) MT=102', color='navy')
plt.xlabel('Energia do nêutron [MeV]')
plt.ylabel('Seção de choque [barns]')
plt.title('Seção de choque Fe-56 (n,xp) (MT=102)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
