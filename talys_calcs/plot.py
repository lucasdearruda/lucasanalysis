import pandas as pd
import matplotlib.pyplot as plt

# Abundâncias relativas
abundancias = {
    "Xs54": 0.058450,
    "Xs56": 0.917540,
    "Xs57": 0.021190,
    "Xs58": 0.002820
}

# Ler o arquivo ignorando linhas que não começam com números
def is_data_line(line):
    return line.strip() and line.strip()[0].isdigit()

with open("testFe0.txt", "r") as f:
    data_lines = [line for line in f if is_data_line(line)]

# Criar DataFrame
from io import StringIO
df = pd.read_csv(StringIO("".join(data_lines)), sep=r'\s+', header=None)

# Verificar número de colunas
print("Número de colunas detectadas:", df.shape[1])

# Assumindo que cada isotopo tem duas colunas: Energia e XS
# Mapear colunas de energia e seção de choque
E = df.iloc[:,0]  # primeira coluna de energia
xs_total = (df.iloc[:,1] * abundancias["Xs54"] +
            df.iloc[:,3] * abundancias["Xs56"] +
            df.iloc[:,5] * abundancias["Xs57"] +
            df.iloc[:,7] * abundancias["Xs58"])

# Plot
plt.figure(figsize=(10,6))
plt.plot(E, xs_total, label="Fe Total Ponderado")
#plt.yscale("lin")
plt.xlabel("Energia (E)")
plt.ylabel("Seção de Choque Ponderada")
plt.title("Seção de Choque Total de Ferro (Ponderada por Abundância)")
plt.grid(True, which="both", ls="--", lw=0.5)
plt.legend()
plt.show()

