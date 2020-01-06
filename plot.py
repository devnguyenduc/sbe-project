import pandas as pd
import matplotlib.pyplot as plt
pn = pd.read_csv("output/detuning_30_Pt.csv")
fn = pd.read_csv("output/detuning_30_Nt.csv")

plt.plot(pn)
plt.show()
plt.plot(fn)
plt.show()