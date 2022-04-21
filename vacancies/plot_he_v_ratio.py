import numpy as np
import matplotlib.pyplot as plt

data_ib = np.genfromtxt("r-ib.txt", delimiter=" ")[::-1]
data_m = np.genfromtxt("r-m.txt", delimiter=" ")[::-1]

indexes = np.where(data_ib[:, 0] < 500e-9)

fig, axs = plt.subplots(3, 1, sharex=True)

axs[0].plot(data_ib[:, 0][indexes], data_ib[:, 1][indexes])
axs[0].set_ylabel(r"$\langle i_b \rangle$ (He)")
axs[0].set_ylim(bottom=0)

axs[1].plot(data_m[:, 0][indexes], data_m[:, 1][indexes])
axs[1].set_ylabel(r"$\langle m_b \rangle$ (V)")
axs[1].set_ylim(bottom=0)


axs[2].plot(data_m[:, 0][indexes], (data_ib[:, 1] / data_m[:, 1])[indexes])
axs[2].set_ylabel(r"$\langle i_b \rangle / \langle m_b \rangle$ (He/V)")
axs[2].set_ylim(bottom=0, top=10)

axs[2].set_xlabel("Depth (m)")

plt.tight_layout()
plt.show()
