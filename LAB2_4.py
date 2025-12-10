import numpy as np

def ladder_analysis(Vs=9.0, Rseries=1_000.0, Rshunt=1_000.0, N=4):

    G = np.zeros((N, N), dtype=float)
    b = np.zeros(N, dtype=float)

    for k in range(1, N + 1):
        i = k - 1

        G[i, i] += 1.0 / Rshunt

        if k == 1:
            G[i, i] += 1.0 / Rseries
            b[i] += Vs / Rseries
        else:
            G[i, i] += 1.0 / Rseries
            G[i, i - 1] -= 1.0 / Rseries

        if k < N:
            G[i, i] += 1.0 / Rseries
            G[i, i + 1] -= 1.0 / Rseries

    V_nodes = np.linalg.solve(G, b)

    I_series = np.zeros(N, dtype=float)
    I_series[0] = (Vs - V_nodes[0]) / Rseries
    for k in range(1, N):
        I_series[k] = (V_nodes[k - 1] - V_nodes[k]) / Rseries

    I_shunt = V_nodes / Rshunt
    Req = Vs / I_series[0]

    return Req, V_nodes, I_series, I_shunt


if __name__ == "__main__":
    Vs = 9.0
    R = 1_000.0
    N = 4

    Req, V_nodes, I_series, I_shunt = ladder_analysis(Vs, R, R, N)

    print(f"Equivalent resistance R_eq = {Req:.2f} ohms")
    print("Node voltages (V1..VN) [V]:")
    print(V_nodes)

    print("Series currents (from left to right) [A]:")
    print(I_series)

    print("Shunt currents (down to ground) [A]:")
    print(I_shunt)
