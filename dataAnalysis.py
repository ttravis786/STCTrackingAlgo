import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_circle(x, y, radius):
    circle = plt.Circle((x, y), radius, fill=False, edgecolor='b')  # b for blue
    plt.gca().add_patch(circle)

df_params = pd.read_csv(r"C++Code/firParams.csv")

print(f"% complet  {100*(df_params.np.sum())/ len(df_params)}")

for i in range(1,101):
    df = pd.read_csv(rf'C++Code/track_{i}pairData.csv')
    df_r = pd.read_csv(rf'C++Code/track_{i}rawData.csv')
    df_p = pd.read_csv(rf'C++Code/track_{i}processedData.csv')
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.plot(df.x, df.y, 'x', label = 'pairPoints')
    #ax.plot(df.x.iloc[0:2], df.y.iloc[0:2], 'x', label = '1')
    #ax.plot(df.x.iloc[2:4], df.y.iloc[2:4], 'x', label = '2')
    #ax.plot(df.x.iloc[4:6], df.y.iloc[4:6], 'x', label = '3')
    ax.plot(df_r.x, df_r.y, label = 'rawPoints')
    ax.plot(df_p.x, df_p.y, '.',  label = 'processedPoints')
    if df_params.iloc[i-1].np == 1:
        ax.set_title('TrackFound')
        x = np.linspace(min(df.x), max(df.x), 10)
        ax.plot(x, df_params.iloc[i-1].m * x + df_params.iloc[i-1].c)
        for x,y,t in zip(df_r.x, df_r.y, df_r.t):
            plot_circle(x, y, t * df_params.iloc[i-1].v/10000)
            print(df_params.iloc[i-1].v/10000)
    else:
        ax.set_title('TrackNotFound')
        for x,y,t in zip(df_r.x, df_r.y, df_r.t):
            plot_circle(x, y, t/500)
    ax.set_aspect('equal', adjustable='box')
    ax.grid()
    fig.tight_layout()
    ax.legend()
    plt.show()

# Keep the script alive while you interact with the plot (if needed)
input("Press [Enter] to exit...")
