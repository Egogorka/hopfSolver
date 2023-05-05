import matplotlib.pyplot as plt
import numpy as np

time = 1

fig, ax = plt.subplots(figsize=(15,15))

def on_press(event):
    global time
    global data
    if event.key == "e":
        time = time + 50
    if event.key == "w":
        time = time - 50
    ax.clear()
    ax.set_xlim(-1, 2+1000)
    ax.set_ylim(0.5, 2.5)
    ax.plot(data[time])
    ax.set_title(f'Time t= {time/5000}')
    # [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % 100 != 0]
    fig.canvas.draw()

fig.canvas.mpl_connect('key_press_event', on_press)

if __name__ == '__main__':
    print("Enter path to data file:")
    path = input()
    # path = "data/cos/data1.csv"

    f = open(path)
    x = f.readline().split(',')
    f.close()
    x.pop(0)

    t = np.loadtxt(path, delimiter=',', skiprows=1, usecols=(1,))
    data = np.loadtxt(path, delimiter=',', skiprows=1, usecols=range(1,1001))

    # print(data[2])
    plt.show()