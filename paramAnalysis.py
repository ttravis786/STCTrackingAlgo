import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv(r'C++Code/firParams.csv')
print(len(df[df.np ==0]))
print(len(df[df.np !=0]))
print(df.columns)
df = df.dropna()
df = df[df.np != 0] 

#df_filtered = df[(df['beta'] - df['beta'].mean()).abs() <= 3 * df['beta'].std()]
plt.hist(df.beta * 180 / np.pi, label='beta',bins=300)
plt.title(f"{np.mean(df.beta * 180 / np.pi)}")
plt.legend()
plt.show()

# also exlucde outliers
v_range = 5
df = df[df.v > df.v.mean() -  v_range]
df = df[df.v <  df.v.mean() +  v_range]
print(len(df))
plt.hist(df.v, label='v', bins=150)
plt.legend()
plt.title(f"{np.mean(df.v)}")
plt.show()
input("Press [Enter] to exit...")

