import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Leer CSV generados por C
# =========================
quat = pd.read_csv("quaternion.csv")
pos = pd.read_csv("position.csv")
vel = pd.read_csv("velocity.csv")
ang = pd.read_csv("ang_velocity.csv")
acc = pd.read_csv("ctrl_accel.csv")
torq = pd.read_csv("ctrl_torque.csv")

# =========================
# Figura 1: Cuaterniones
# =========================
plt.figure(1)
plt.plot(quat["time"], quat["q0"], label="q0")
plt.plot(quat["time"], quat["q1"], label="q1")
plt.plot(quat["time"], quat["q2"], label="q2")
plt.plot(quat["time"], quat["q3"], label="q3")

qnorm = np.sqrt(
    quat["q0"]**2 + quat["q1"]**2 + quat["q2"]**2 + quat["q3"]**2
)
plt.plot(quat["time"], qnorm, label="||q||")

plt.ylabel("Quaternion sequence")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

# =========================
# Figura 2: Posición
# =========================
plt.figure(2)
plt.plot(pos["time"], pos["x"], label="x")
plt.plot(pos["time"], pos["y"], label="y")
plt.plot(pos["time"], pos["z"], label="z")

plt.ylabel("Position sequence [m]")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

# =========================
# Figura 3: Velocidad
# =========================
plt.figure(3)
plt.plot(vel["time"], vel["vx"], label="vx")
plt.plot(vel["time"], vel["vy"], label="vy")
plt.plot(vel["time"], vel["vz"], label="vz")

plt.ylabel("Velocity sequence [m/s]")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

# =========================
# Figura 4: Velocidad angular
# =========================
plt.figure(4)
plt.plot(ang["time"], ang["wx"], label="wx")
plt.plot(ang["time"], ang["wy"], label="wy")
plt.plot(ang["time"], ang["wz"], label="wz")

plt.ylabel("Angular velocity sequence [rad/s]")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

# =========================
# Figura 5: Aceleraciones de control
# =========================
plt.figure(5)
plt.step(acc["time"], acc["ax"], label="ax", where="post")
plt.step(acc["time"], acc["ay"], label="ay", where="post")
plt.step(acc["time"], acc["az"], label="az", where="post")

plt.ylabel("Control acceleration sequence [m/s^2]")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

# =========================
# Figura 6: Torques de control
# =========================
plt.figure(6)
plt.step(torq["time"], torq["taux"], label="taux", where="post")
plt.step(torq["time"], torq["tauy"], label="tauy", where="post")
plt.step(torq["time"], torq["tauz"], label="tauz", where="post")

plt.ylabel("Control torque sequence [Nm]")
plt.xlabel("Time [s]")
plt.gca().tick_params(labelsize=10)
plt.legend()
plt.grid(True)

plt.show()